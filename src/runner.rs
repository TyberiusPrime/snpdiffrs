use crate::chunked_genome;
use crate::chunked_genome::Chunk;
use crate::consts::*;
use crate::coverage::{Coverage, ResultRow};

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::Mutex;

use itertools::{iproduct, Itertools};
use progressing::{mapping::Bar as MappingBar, Baring};
use rust_htslib::bam;
use serde::Deserialize;

#[allow(unused_imports)]
use slog::debug;

use crossbeam::thread;

fn default_quality_threshold() -> u8 {
    15u8
}

fn default_blocksize() -> usize {
    50_000_000
}

#[derive(Deserialize, Debug)]
pub struct RunConfig {
    output_dir: String,
    chromosomes: Option<Vec<String>>,
    samples: HashMap<String, Vec<String>>,
    #[serde(default = "default_blocksize")]
    block_size: usize,
    ncores: Option<u32>,
    #[serde(default = "default_quality_threshold")]
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: Option<f32>,
}

impl RunConfig {
    fn samples_and_input_filenames(&self) -> Vec<(String, Vec<PathBuf>)> {
        let mut temp: Vec<(&String, &Vec<String>)> = self.samples.iter().collect();
        temp.sort();
        temp.iter()
            .map(|(name, bam_files)| (name.to_string(), bam_files.iter().map(|x| x.into()).collect::<Vec<PathBuf>>()))
            .collect()
    }
}

pub fn snp_diff_from_toml(input: &str) -> Result<(), ()> {
    let config: RunConfig = toml::from_str(input).unwrap();
    run_snp_diff(config)
}

pub fn run_snp_diff(config: RunConfig) -> Result<(), ()> {
    for filenames in config.samples.values() {
        for filename in filenames {
            if !Path::new(filename).exists() {
                panic!("File did not exist {}", filename);
            }
        }
    }
    let output_dir = Path::new(&config.output_dir);
    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).unwrap();
    }

    let first_bam = config
        .samples
        .values()
        .next()
        .expect("No samples")
        .iter()
        .next()
        .unwrap();
    let first_bam = bam::IndexedReader::from_path(first_bam).unwrap();

    let chunks = chunked_genome::ChunkedGenome::new(first_bam, &config.chromosomes);
    let block_size = config.block_size;
    let runner = NtoNRunner::new(config, chunks.iter(block_size).collect());
    runner.run();
    //files are now closed. rename them all
    Ok(())
}

fn write_results(
    output_handle: &mut std::io::BufWriter<std::fs::File>,
    chr: &str,
    chunk_start: u32,
    rows: Vec<ResultRow>,
) {
    if output_handle.seek(SeekFrom::Current(0)).unwrap() == 0 {
        //nothing written yet
        writeln!(
            output_handle,
            "\tchr\tpos\tscore\tA_A\tB_A\tA_C\tB_C\tA_G\tB_G\tA_T\tB_T\thaplotypeA\thaplotypeB"
        )
        .unwrap();
    }
    for (number, row) in rows.iter().enumerate() {
        writeln!(
            output_handle,
            "{}\t{}\t{}\t{:.13}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            number,
            chr,
            row.relative_pos + chunk_start,
            row.score,
            row.count_self_a,
            row.count_other_a,
            row.count_self_c,
            row.count_other_c,
            row.count_self_g,
            row.count_other_g,
            row.count_self_t,
            row.count_other_t,
            haplotype_const_to_str(row.haplotype_self),
            haplotype_const_to_str(row.haplotype_other),
        )
        .unwrap();
    }
}

fn haplotype_const_to_str(haplotype: u8) -> &'static str {
    match haplotype as usize {
        AA => "AA",
        AC => "AC",
        AG => "AG",
        AT => "AT",
        CC => "CC",
        CG => "CG",
        CT => "CT",
        GG => "GG",
        GT => "GT",
        TT => "TT",
        NN => "NN",
        _ => "??",
    }
}

#[allow(dead_code)]
fn get_logger() -> slog::Logger {
    use slog::{o, Drain};

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    slog::Logger::root(drain, o!())
}

#[derive(Debug)]
struct PayloadTodoLoadCoverage{sample_id: usize,
chunk: Chunk,
chunk_id: usize
}

enum MsgTodo {
    LoadCoverage(PayloadTodoLoadCoverage),
    CalcSnps(usize, [usize; 2], Arc<Coverage>, Arc<Coverage>, Chunk),
    OutputResult(Vec<ResultRow>, Chunk, Arc<Mutex<BufWriter<File>>>),
    Quit,
}

impl std::fmt::Debug for MsgTodo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MsgTodo::LoadCoverage(_payload) => {
                f.write_fmt(format_args!("MsgTodo::LoadCoverage"))
            }
            MsgTodo::CalcSnps(_, _, _, _, _) => f.write_fmt(format_args!("MsgTodo:Calc_Snps")),
            MsgTodo::OutputResult(_, _, _) => f.write_fmt(format_args!("MsgTodo:OutputResult")),
            MsgTodo::Quit => f.write_fmt(format_args!("MsgTodo:Quit")),
        }
    }
}

struct PayloadDoneLoadCoverage {
    chunk_id: usize,
    sample_id: usize,
    coverage: Coverage,
    chunk: Chunk,
}

enum MsgDone {
    LoadedCoverage(PayloadDoneLoadCoverage), //chunk, sample_id,
    CalculatedSnps(usize, [usize; 2], Vec<ResultRow>, Chunk),
    OutputDone,
    QuitDone,
}

impl std::fmt::Debug for MsgDone {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MsgDone::LoadedCoverage(payload) => {
                f.write_fmt(format_args!(
                    "MsgDone::LoadedCoverage(chunk_id: {}, sample_id: {})",
                    payload.chunk_id, payload.sample_id
                ))
            }
            MsgDone::CalculatedSnps(chunk_id, pair, result_rows, _chunk) => {
                f.write_fmt(format_args!(
                    "MsgDone::CalculatedSnp(c: {}, pair: {:?}, len(results)=={}",
                    chunk_id,
                    pair,
                    result_rows.len()
                ))
            }
            MsgDone::OutputDone => f.write_fmt(format_args!("JobResult::OutputDone")),
            MsgDone::QuitDone => f.write_fmt(format_args!("JobResult::QuitDone")),
        }
    }
}

struct NtoNRunner {
    config: RunConfig,
    chunks: Vec<Chunk>,
}

#[derive(Debug)]
struct Block {
    chunk_id: usize,
    sample_id: usize,
    coverage: Arc<Coverage>,
    used: usize,
}

impl NtoNRunner {
    fn new(config: RunConfig, chunks: Vec<Chunk>) -> NtoNRunner {
        NtoNRunner { config, chunks }
    }

    fn run(self) {
        let log = get_logger();

        //input
        let samples_and_input_filenames= self.config.samples_and_input_filenames();
        let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) = samples_and_input_filenames.into_iter().unzip();
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_filenames.into_iter().map(Arc::new).collect();
        let block_iterator: Vec<_> =
            iproduct!(self.chunks.iter().enumerate(), 0..input_filenames.len())
                .enumerate()
                .collect();
        let mut block_iterator: Vec<_> = block_iterator.into_iter().rev().collect();
        let quality_threshold = self.config.quality_threshold;
        let filter_homo_polymer_threshold = self.config.filter_homo_polymer_threshold;
        let min_score = self.config.min_score.unwrap_or(50.0);

        //output
        let output_dir = Path::new(&self.config.output_dir);
        let pairs: Vec<[usize; 2]> = (0..input_filenames.len())
            .combinations(2)
            .map(|x| [x[0], x[1]])
            .collect();
        let output_filenames: Vec<_> = pairs
            .iter()
            .map(|x| output_dir.join(format!("{}_vs_{}.tsv.tmp",
                                             &samples[x[0]],
                                             &samples[x[1]],
                                             )))
            .collect();
        let outputs: Vec<_> = output_filenames
            .iter()
            .map(|filename| std::fs::File::create(filename).unwrap())
            .map(std::io::BufWriter::new)
            .map(|buf| Arc::new(Mutex::new(buf)))
            .collect();

        let outputs: HashMap<[usize; 2], _> = pairs.iter().copied().zip(outputs).collect();

        //tunables
        let ncores = self.config.ncores.unwrap_or(num_cpus::get() as u32);
        let max_concurrent_blocks = input_filenames.len() + ncores as usize; //we need at least one per input_filename to create all pairs. And some bonus is helpful to avoid delays.
         //debug!(log, "ncores: {}", ncores);

        let mut load_progress = MappingBar::with_range(
            0,
            block_iterator.len() + (input_filenames.len().pow(2) / 2 * self.chunks.len()),
        )
        .timed();
        load_progress.set_len(40);
        load_progress.set(0usize);
        print!("{}", load_progress);

        let (todo_sender, todo_receiver) = crossbeam::crossbeam_channel::unbounded();
        for ii in 0..max_concurrent_blocks {
            if block_iterator.is_empty() {
                if ii == 0 {
                    panic!("Nothing to do?");
                } else {
                    break;
                }
            }
            let (_, ((chunk_id, chunk), sample_id)) = block_iterator.pop().unwrap();
            todo_sender
                .send(MsgTodo::LoadCoverage(
                        PayloadTodoLoadCoverage{
                            sample_id: sample_id,
                            chunk: chunk.clone(),
                            chunk_id: chunk_id

                        }
                ))
                .expect("Could not send initial loads");
        }

        //runtime stuff
        let (result_sender, result_receiver) = crossbeam::crossbeam_channel::unbounded();
        let block_iterator = Arc::new(Mutex::new(block_iterator));
        let expected_number_of_uses_per_block = { input_filenames.len() -1};  //in each pair we count both blocks, but we don't compare a block to itself.
        debug!(log, "Expected number of uses for each block: {}", expected_number_of_uses_per_block);
        let mut blocks: Vec<Block> = Vec::new(); //the currently loaded blocks.

        let mut quit_counter = 0;
        thread::scope(|s| {
        let block_iterator = block_iterator.clone();
        for nn in 0..ncores {
            let thread_recv = todo_receiver.clone();
            let result_sender = result_sender.clone();
            let log = log.clone();
            let input_filenames = input_filenames.clone(); // todo: replace with arc,  I suppose.
            s.spawn(move |_s2| {
                //so this is what happens in the worker threads.
                for el in thread_recv {
                    match el {
                        MsgTodo::LoadCoverage(payload) => {
                            //debug!(log, "exc: Load_coverage: c:{} i:{}", chunk_no, input_no);
                            debug!(log, "Loading {} {}", payload.chunk.chr, payload.chunk.start);
                            result_sender
                                .send(MsgDone::LoadedCoverage(
                                        PayloadDoneLoadCoverage{
                                            chunk_id: payload.chunk_id,
                                            sample_id: payload.sample_id,
                                            coverage:
                                            Coverage::from_bams(&input_filenames[payload.sample_id].as_ref().iter().map(|x| x.as_path()).collect(),
                                                        payload.chunk.tid,
                                                        payload.chunk.start,
                                                        payload.chunk.stop,
                                                        quality_threshold,
                                                        &filter_homo_polymer_threshold,
                                                        ),
                                            chunk: payload.chunk,
                                        }
                                ))
                                .expect("Could not send LoadedCoverage reply");
                        }
                        MsgTodo::CalcSnps(chunk_id, pair, cov, other_cov, chunk) => {
                            let result: Vec<ResultRow> = other_cov.score_differences(cov.as_ref(), min_score, chunk.start);
                            //debug!(log, "exc: calc snps c: {}, pair: {:?}, result_len: {}", chunk_id, pair, result.len());
                            debug!(log, "Calced {}, {}" ,chunk.chr, chunk.start);
                            result_sender.send(
                                MsgDone::CalculatedSnps(chunk_id, pair, result, chunk)).expect("Could not send CalculatedSnps");
                        }
                        MsgTodo::OutputResult(result_rows, chunk, output) => {
                            //debug!(log, "exc: output results for chr: {}, start: {}, count: {}, first entry: {}", chunk.chr, chunk.start, result_rows.len(), if result_rows.is_empty() {0} else {result_rows[0].relative_pos});
                            debug!(log, "writing {}, {}", chunk.chr, chunk.start);
                            write_results(
                                &mut output.lock().unwrap(),
                                &chunk.chr,
                                0, result_rows);
                            result_sender.send(MsgDone::OutputDone).expect("could not signal OutputDone");
                        }
                        MsgTodo::Quit => {
                            debug!(log, "exc: quit: {}", nn);
                            result_sender.send(MsgDone::QuitDone).expect("Could not send QuitDone");
                            return;
                        }
                    }
                }
            });
        }
        s.spawn(move |_s2| {
            for res in result_receiver.iter() {
                debug!(log, "Received a result: {:?}", res);
                match res {
                    MsgDone::LoadedCoverage(payload) => {
                        //Todo: optimize by not putting empty blocks into the pool
                        //debug!(log, "Loadeded Coverage {} {}", chunk_id, sample_id);
                        debug!(log, "loaded {} {}", payload.chunk.chr, payload.chunk.start);
                        load_progress.add(1usize);
                        if load_progress.has_progressed_significantly() {
                            print!("\r{}", load_progress);
                        }

                        let cov = Arc::new(payload.coverage);
                        let used = expected_number_of_uses_per_block;
                        if !blocks.is_empty() {
                            for block in
                                blocks.iter_mut()
                            {
                                //debug!(log, "checking block other_chunk_id {}: other_sample_id, {}", block.chunk_id, block.sample_id);
                                if block.chunk_id == payload.chunk_id{
                                    let (a,b, cov_a, cov_b) =
                                        if block.sample_id < payload.sample_id {
                                            (block.sample_id, payload.sample_id, cov.clone(), block.coverage.clone())
                                        }
                                        else  {
                                            (payload.sample_id, block.sample_id, block.coverage.clone(), cov.clone())
                                        };
                                    debug!( log, "Sending pair {} {} for chunk {} {}", a, b, payload.chunk.chr, payload.chunk.start);
                                    todo_sender
                                        .send(MsgTodo::CalcSnps(
                                                payload.chunk_id,
                                                [a, b],
                                                cov_a,
                                                cov_b,
                                                payload.chunk.clone(),
                                                ))
                                        .expect("Could not send CalcSnps");
                                    //block.used += 1;
                                    //used += 1;
                                }
                            }
                        }
                        blocks.push(Block{chunk_id: payload.chunk_id, sample_id: payload.sample_id, coverage: cov, used});
                        },
                    MsgDone::CalculatedSnps(chunk_id, pair, result_rows, chunk) =>{
                            debug!(log, "Received CalculatedSnps for c:{}, pair: {:?}", chunk_id, pair);
                            load_progress.add(1usize);
                            if load_progress.has_progressed_significantly() {
                                print!("\r{}", load_progress);
                            }
                            todo_sender.send(MsgTodo::OutputResult(result_rows, chunk,
                                                                outputs.get(&pair).expect("Could not find output for pair").clone())
                                    ).expect("Could not send OutputResult");
                            let mut to_delete = Vec::new();
                            let mut found = 0;
                            for (ii, block) in blocks.iter_mut().enumerate() {
                                if (block.chunk_id == chunk_id) && ((block.sample_id == pair[0]) || (block.sample_id == pair[1])) {
                                    debug!(log, "Accepted block: {:?}", block);
                                    block.used -= 1;
                                    found += 1;
                                    if block.used == 0{
                                        debug!(log, "deleting block c: {}, i: {}; ii: {}", block.chunk_id, block.sample_id, ii);
                                        to_delete.push(ii);
                                    }
                                }
                            }
                            if found != 2 {
                                panic!("you made a mistake and we got back a block that no longer exists: c: {}, pair: {:?}, found: {}", chunk_id, pair, found);
                            }
                            for ii in to_delete.iter().rev() {
                                debug!(log, "deleting block at {}", ii);
                                blocks.remove(*ii);
                            }
                        if blocks.len() < max_concurrent_blocks {
                                for _ in blocks.len()..max_concurrent_blocks {
                                    let next = block_iterator.lock().expect("Could not lock block_iterator").pop(); 
                                    match next {
                                        Some((_, ((chunk_no, chunk), sample_id))) => {
                                            todo_sender
                                                .send(MsgTodo::LoadCoverage(
                                                        PayloadTodoLoadCoverage{
                                                            sample_id: sample_id,
                                                            chunk: chunk.clone(),
                                                            chunk_id: chunk_no,
                                                        }
                                                ))
                                                .unwrap();
                                    },
                                    None => {},
                                    };
                                }
                            }
                        },
                        MsgDone::OutputDone => {
                            if block_iterator.lock().unwrap().is_empty() && blocks.is_empty() { //no more work incoming, and everything done.
                                //debug!(log, "block_iterator empty - preparing to leave");
                                for _nn in 0..ncores {
                                    todo_sender.send(MsgTodo::Quit).unwrap();
                                }
                            }


                        }
                        MsgDone::QuitDone => {
                            quit_counter += 1;
                            if quit_counter == ncores {
                                break;
                            }

                        }
                    }
                }
                //everything processed
                drop(outputs);
        });
        })
        .expect("An error occured in one of the scoped threads");

        output_filenames
            .iter()
            .map(|tmp_filename| {
                std::fs::rename(
                    tmp_filename,
                    tmp_filename.with_file_name(tmp_filename.file_stem().unwrap()),
                )
                .unwrap()
            })
            .for_each(|_| {});
    }
}

mod test {

    #[test]
    fn test_sample_data() {
        let toml = "
output_dir = 'tests/test_sample_data'
[samples]
    A = ['sample_data/sample_a.bam']
    B = ['sample_data/sample_b.bam']
";
        use std::path::Path;
        let output_dir = Path::new("tests/test_sample_data");
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        super::snp_diff_from_toml(toml).unwrap();
        let should =
            std::fs::read_to_string("sample_data/marsnpdiff_sample_a_vs_sample_b.tsv").unwrap();
        let should = should.replace(".0", "");
        let actual = std::fs::read_to_string("tests/test_sample_data/A_vs_B.tsv").unwrap();
        assert_eq!(should, actual)
    }
    #[test]
    fn test_sample_n_to_n() {
        let test_path = "tests/test_sample_data_n_to_n";
        let toml = format!("
output_dir = '{}'
[samples]
    A = ['sample_data/sample_a.bam'] # mixed up order is intentional
    C = ['sample_data/sample_a.bam']
    D = ['sample_data/sample_b.bam']
    B = ['sample_data/sample_b.bam']
", test_path);
        use std::path::Path;
        let output_dir = Path::new(test_path);
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        super::snp_diff_from_toml(&toml).unwrap();
        let should =
            std::fs::read_to_string("sample_data/marsnpdiff_sample_a_vs_sample_b.tsv").unwrap();
        let should = should.replace(".0", "");

        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_B.tsv").unwrap();
        assert_eq!(should, actual);
        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_D.tsv").unwrap();
        assert_eq!(should, actual);

        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_C.tsv").unwrap();
        assert_eq!(actual.matches("\n").count(), 1);

        let actual = std::fs::read_to_string(test_path.to_owned() + "/B_vs_D.tsv").unwrap();
        assert_eq!(actual.matches("\n").count(), 1);
        std::fs::remove_dir_all(output_dir).ok();


    }


}
