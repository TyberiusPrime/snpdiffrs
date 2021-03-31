use crate::chunked_genome::{Chunk, ChunkedGenome};
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

use crossbeam::thread;
#[allow(unused_imports)]
use slog::debug;

use super::{all_files_exists, ensure_output_dir, get_logger, haplotype_const_to_str, RunConfig};

pub fn run_n_to_n(config: RunConfig) -> Result<(), ()> {
    all_files_exists(&config.samples);
    ensure_output_dir(&Path::new(&config.output_dir));

    let chunks = ChunkedGenome::new(config.first_bam(), &config.chromosomes);
    let block_size = config.block_size;
    let runner = NtoNRunner::new(config, chunks.iter(block_size).collect());
    runner.run();
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

#[derive(Debug)]
struct PayloadTodoLoadCoverage {
    sample_id: usize,
    chunk: Chunk,
    chunk_id: usize,
}

#[derive(Debug)]
struct PayloadTodoCalcSnps {
    chunk_id: usize,
    sample_pair: [usize; 2],
    coverage_a: Arc<Coverage>,
    coverage_b: Arc<Coverage>,
    chunk: Chunk,
}

#[derive(Debug)]
struct PayloadTodoOutputResult {
    result_rows: Vec<ResultRow>,
    chunk: Chunk,
    output: Arc<Mutex<BufWriter<File>>>,
}

#[derive(Debug)]
enum MsgTodo {
    LoadCoverage(PayloadTodoLoadCoverage),
    CalcSnps(PayloadTodoCalcSnps),
    OutputResult(PayloadTodoOutputResult),
}

#[derive(Debug)]
struct PayloadDoneLoadCoverage {
    chunk_id: usize,
    sample_id: usize,
    coverage: Coverage,
    chunk: Chunk,
}

#[derive(Debug)]
struct PayloadDoneCalcSnps {
    chunk_id: usize,
    sample_pair: [usize; 2],
    result_rows: Vec<ResultRow>,
    chunk: Chunk,
}

#[derive(Debug)]
enum MsgDone {
    LoadCoverage(PayloadDoneLoadCoverage), //chunk, sample_id,
    CalcSnps(PayloadDoneCalcSnps),
    OutputDone,
}

#[allow(clippy::type_complexity)]
struct NtoNRunner {
    //config: RunConfig,
    //chunks: Vec<Chunk>,
    samples: Vec<String>,
    input_filenames: Vec<Arc<Vec<PathBuf>>>,
    output_filenames: Vec<PathBuf>,
    outputs: Option<HashMap<[usize; 2], Arc<Mutex<std::io::BufWriter<std::fs::File>>>>>,
    preprocessed_dir: Option<PathBuf>,

    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: f32,

    ncores: u32,
    max_concurrent_blocks: usize,

    chunks: Vec<Chunk>,
    block_size: usize,
}

#[derive(Debug)]
struct Block {
    chunk_id: usize,
    sample_id: usize,
    coverage: Arc<Coverage>,
    remaining_uses: usize,
}

impl NtoNRunner {
    fn new(config: RunConfig, chunks: Vec<Chunk>) -> NtoNRunner {
        //input
        let samples_and_input_filenames = config.samples_and_input_filenames();
        let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) =
            samples_and_input_filenames.into_iter().unzip();
        if samples.len() < 2 {
            panic!("Not enough samples supplied");
        }
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_filenames.into_iter().map(Arc::new).collect();

        //output
        let output_dir = Path::new(&config.output_dir);
        let pairs: Vec<[usize; 2]> = (0..input_filenames.len())
            .combinations(2)
            .map(|x| [x[0], x[1]])
            .collect();
        let output_filenames: Vec<_> = pairs
            .iter()
            .map(|x| output_dir.join(format!("{}_vs_{}.tsv.tmp", &samples[x[0]], &samples[x[1]],)))
            .collect();
        let outputs: Vec<_> = output_filenames
            .iter()
            .map(|filename| std::fs::File::create(filename).unwrap())
            .map(std::io::BufWriter::new)
            .map(|buf| Arc::new(Mutex::new(buf)))
            .collect();
        let outputs: HashMap<[usize; 2], _> = pairs.iter().copied().zip(outputs).collect();

        let ncores = config.ncores.unwrap_or(num_cpus::get() as u32);
        let max_concurrent_blocks =
            (input_filenames.len() + ncores as usize).max(ncores as usize * 3); //we need at least one per input_filename to create all pairs. And some bonus is helpful to avoid delays.

        NtoNRunner {
            samples,
            input_filenames,
            output_filenames,
            outputs: Some(outputs),
            quality_threshold: config.quality_threshold,
            filter_homo_polymer_threshold: config.filter_homo_polymer_threshold,
            min_score: config.min_score.unwrap_or(50.0),
            ncores,
            max_concurrent_blocks,
            chunks,
            preprocessed_dir: config.preprocessed_dir.map(|x| PathBuf::from(&x)),
            block_size: config.block_size,
        }
    }

    #[allow(clippy::type_complexity)]
    fn _push_next_block(
        block_iterator: &mut Vec<(usize, ((usize, Chunk), usize))>,
        todo_sender: &crossbeam::Sender<MsgTodo>,
    ) {
        if let Some((_, ((chunk_id, chunk), sample_id))) = block_iterator.pop() {
            todo_sender
                .send(MsgTodo::LoadCoverage(PayloadTodoLoadCoverage {
                    sample_id,
                    chunk,
                    chunk_id,
                }))
                .expect("Could not send initial loads");
        }
    }

    fn run(self) {
        let log = get_logger();

        let chunk_count = self.chunks.len();
        let block_iterator: Vec<_> = iproduct!(
            self.chunks.clone().into_iter().enumerate(),
            0..self.input_filenames.len()
        )
        .enumerate()
        .collect();
        let mut block_iterator: Vec<_> = block_iterator.into_iter().rev().collect();
        //tunables
        //debug!(log, "ncores: {}", ncores);

        let expected_number_of_uses_per_block = { self.input_filenames.len() - 1 }; //in each pair we count both blocks, but we don't compare a block to itself.
        let expected_number_of_outputs =
            self.input_filenames.len() * (self.input_filenames.len() - 1) / 2 * chunk_count;

        let mut load_progress = MappingBar::with_range(
            0,
            block_iterator.len()
                + block_iterator.len() * expected_number_of_uses_per_block * chunk_count,
        )
        .timed();
        load_progress.set_len(40);
        load_progress.set(0usize);
        print!("{}", load_progress);
        if block_iterator.is_empty() {
            panic!("Nothing to do?");
        }

        let (todo_sender, todo_receiver) = crossbeam::crossbeam_channel::unbounded();
        let mut _self = Arc::new(self);
        for _ in 0.._self.max_concurrent_blocks {
            Self::_push_next_block(&mut block_iterator, &todo_sender);
        }

        //runtime stuff
        let (result_sender, result_receiver) = crossbeam::crossbeam_channel::unbounded();
        let block_iterator = Arc::new(Mutex::new(block_iterator));
        debug!(
            log,
            "Expected number of uses for each block: {}", expected_number_of_uses_per_block
        );
        debug!(
            log,
            "Expected number of output blocks: {}", expected_number_of_outputs
        );
        let mut blocks: Vec<Block> = Vec::new(); //the currently loaded blocks.
        let remaining_outputs = Mutex::new(expected_number_of_outputs);

        thread::scope(|s| {
        let block_iterator = block_iterator.clone();
        let _self = _self.clone();
        for nn in 0.._self.ncores {
            let thread_recv = todo_receiver.clone();
            let result_sender = result_sender.clone();
            let log = log.clone();
            let _self = _self.clone();
            s.spawn(move |_s2| {
                //so this is what happens in the worker threads.
                for el in thread_recv {
                    match el {
                        MsgTodo::LoadCoverage(payload) => {
                            //debug!(log, "exc: Load_coverage: c:{} i:{}", chunk_no, input_no);
                            debug!(log, "Loading {} {}", payload.chunk.chr, payload.chunk.start);
                            let cov: Option<Coverage> = match _self.preprocessed_dir.as_ref() {
                                Some(pd) => {
                                    let full_pd = pd.join(format!(
                                        "{}_bs={}_q={}_homo_polymer={}.snpdiff_block",_self.samples[payload.sample_id],
                                        _self.block_size,
                                        _self.quality_threshold,
                                        match _self.filter_homo_polymer_threshold {
                                                                Some(x) => format!("{}", x),
                                                                None => "None".to_string(),
                                    }));

                                    let filename = full_pd.join(
                                        format!("{}_{}.snpdiffrs_block",
                                            payload.chunk.chr, payload.chunk.start));
                                    if !filename.exists() {
                                        println!("Could not load preprocessed data: {:?}", filename);
                                        None
                                    } else {
                                        Coverage::from_preprocessed(&filename)
                                    }
                                },
                                None => None
                            };
                            let cov: Coverage = match cov {
                                Some(cov) => cov,
                                None => Coverage::from_bams(&_self.input_filenames[payload.sample_id].as_ref().iter().map(|x| x.as_path()).collect::<Vec<_>>(),
                                        payload.chunk.tid,
                                        payload.chunk.start,
                                        payload.chunk.stop,
                                        _self.quality_threshold,
                                        &_self.filter_homo_polymer_threshold,
                                )
                            };


                            result_sender
                                .send(MsgDone::LoadCoverage(
                                        PayloadDoneLoadCoverage{
                                            chunk_id: payload.chunk_id,
                                            sample_id: payload.sample_id,
                                            coverage: cov,
                                            chunk: payload.chunk,
                                        }
                                ))
                                .expect("Could not send LoadCoverage reply");
                        },

                        MsgTodo::CalcSnps(payload) => {
                            let result: Vec<ResultRow> = payload.coverage_b.score_differences(
                                payload.coverage_a.as_ref(),
                                _self.min_score,
                                payload.chunk.start);
                            //let result: Vec<ResultRow> = Vec::new();
                            //debug!(log, "exc: calc snps c: {}, pair: {:?}, result_len: {}", chunk_id, pair, result.len());
                            debug!(log, "Calced {}, {}" ,payload.chunk.chr, payload.chunk.start);
                            result_sender.send(
                                MsgDone::CalcSnps(PayloadDoneCalcSnps{
                                    chunk_id: payload.chunk_id,
                                    sample_pair: payload.sample_pair,
                                    result_rows: result,
                                    chunk: payload.chunk})).expect("Could not send CalcSnps");
                        }

                        MsgTodo::OutputResult(payload) => {
                            //debug!(log, "exc: output results for chr: {}, start: {}, count: {}, first entry: {}", chunk.chr, chunk.start, result_rows.len(), if result_rows.is_empty() {0} else {result_rows[0].relative_pos});
                            debug!(log, "writing {}, {} {}", payload.chunk.chr, payload.chunk.start, payload.result_rows.len());
                            write_results(
                                &mut payload.output.lock().unwrap(),
                                &payload.chunk.chr,
                                0, payload.result_rows);
                            result_sender.send(MsgDone::OutputDone).expect("could not signal OutputDone");
                        }
                    }
                }
            });
        }
        s.spawn(move |_s2| {
            //and this is the managing thread
            //needs to be a separate spawn so that the panics! work as intended
            for res in result_receiver.iter() {
                debug!(log, "Received a result: {:?}", res);
                match res {
                    MsgDone::LoadCoverage(payload) => {
                        //Todo: optimize by not putting empty blocks into the pool
                        //debug!(log, "Loadeded Coverage {} {}", chunk_id, sample_id);
                        debug!(log, "loaded {} {} {}", payload.chunk.chr, payload.chunk.start, payload.coverage.sum());
                        load_progress.add(1usize);
                        if load_progress.has_progressed_significantly() {
                            print!("\r{}", load_progress);
                        }

                        let cov = Arc::new(payload.coverage);
                        let remaining_uses = expected_number_of_uses_per_block;
                        if !blocks.is_empty() {
                            for block in
                                blocks.iter_mut()
                            {
                                if block.chunk_id == payload.chunk_id{
                                    let (a,b, cov_a, cov_b) =
                                        if block.sample_id < payload.sample_id {
                                            (block.sample_id, payload.sample_id, cov.clone(), block.coverage.clone())
                                        }
                                        else  {
                                            (payload.sample_id, block.sample_id, block.coverage.clone(), cov.clone())
                                        };
                                    debug!( log, "Sending pair {} {} for chunk {} {}", a, b, payload.chunk.chr, payload.chunk.start);
                                    //panic!("how much is the memory");
                                    todo_sender
                                        .send(MsgTodo::CalcSnps(
                                                PayloadTodoCalcSnps{
                                                    chunk_id:payload.chunk_id,
                                                    sample_pair: [a, b],
                                                    coverage_a: cov_a,
                                                    coverage_b: cov_b,
                                                    chunk: payload.chunk.clone()
                                                }
                                                ))
                                        .expect("Could not send CalcSnps");
                                }
                            }
                        }
                        blocks.push(Block{chunk_id: payload.chunk_id, sample_id: payload.sample_id, coverage: cov, remaining_uses});
                    },

                    MsgDone::CalcSnps(payload) => {
                            debug!(log, "Received CalcSnps for c:{}, pair: {:?}", payload.chunk_id, payload.sample_pair);
                            load_progress.add(1usize);
                            if load_progress.has_progressed_significantly() {
                                print!("\r{}", load_progress);
                            }
                            todo_sender.send(MsgTodo::OutputResult(
                                    PayloadTodoOutputResult{
                                        result_rows: payload.result_rows,
                                        chunk: payload.chunk,
                                        output: _self.outputs.as_ref().unwrap().get(&payload.sample_pair).expect("Could not find output for pair").clone()
                                    })).expect("Could not send OutputResult");
                            let mut to_delete = Vec::new();
                            let mut found = 0;
                            for (ii, block) in blocks.iter_mut().enumerate() {
                                if (block.chunk_id == payload.chunk_id) && ((block.sample_id == payload.sample_pair[0]) || (block.sample_id == payload.sample_pair[1])) {
                                    debug!(log, "Accepted block: {:?}", block);
                                    block.remaining_uses -= 1;
                                    found += 1;
                                    if block.remaining_uses == 0{
                                        debug!(log, "deleting block c: {}, i: {}; ii: {}", block.chunk_id, block.sample_id, ii);
                                        to_delete.push(ii);
                                    }
                                }
                            }
                            if found != 2 {
                                panic!("a mistake happend and we got back a block that no longer exists: c: {}, pair: {:?}, found: {}", payload.chunk_id, payload.sample_pair, found);
                            }
                            for ii in to_delete.iter().rev() {
                                debug!(log, "deleting block at {} - arc count strong: {} weak; {}", ii, 
                                    Arc::strong_count(&blocks[*ii].coverage),
                                    Arc::weak_count(&blocks[*ii].coverage),
                                );
                                blocks.remove(*ii);
                            }
                        if blocks.len() < _self.max_concurrent_blocks {
                                for _ in blocks.len().._self.max_concurrent_blocks {
                                Self::_push_next_block(&mut block_iterator.lock().unwrap(), &todo_sender);
                                }
                            }
                        },
                    MsgDone::OutputDone => {
                        let mut rr = remaining_outputs.lock().unwrap();
                            *rr -= 1;
                                debug!(log, "MsgDone:OutputDone - remaining outputs {}", *rr);
                            if *rr == 0 {
                                debug!(log, "all expected outputs done");
                                break
                            };
                    }
                }
            }
                //everything processed
        });
        })
        .expect("An error occured in one of the scoped threads");
        let log = get_logger();
        debug!(log, "Left loop");

        drop(Arc::get_mut(&mut _self).unwrap().outputs.take());
        _self
            .output_filenames
            .iter()
            .map(|tmp_filename| {
                std::fs::rename(
                    tmp_filename,
                    tmp_filename.with_file_name(tmp_filename.file_stem().unwrap()),
                )
                .unwrap()
            })
            .for_each(|_| {});
        debug!(log, "Done with dropping");
    }
}

mod test {
    #[cfg(test)]
    use super::super::snp_diff_from_toml;

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
        snp_diff_from_toml(toml).unwrap();
        let should =
            std::fs::read_to_string("sample_data/marsnpdiff_sample_a_vs_sample_b.tsv").unwrap();
        let should = should.replace(".0", "");
        let actual = std::fs::read_to_string("tests/test_sample_data/A_vs_B.tsv").unwrap();
        assert_eq!(should, actual)
    }
    #[test]
    #[should_panic]
    fn test_same_data() {
        let toml = "
output_dir = 'tests/test_same_data'
[samples]
    A = ['sample_data/sample_a.bam']
    A = ['sample_data/sample_a.bam']
";
        use std::path::Path;
        let output_dir = Path::new("tests/test_same_data");
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        snp_diff_from_toml(toml).unwrap();
    }

    #[test]
    fn test_sample_n_to_n() {
        let test_path = "tests/test_sample_data_n_to_n";
        let toml = format!(
            "
output_dir = '{}'
[samples]
    A = ['sample_data/sample_a.bam'] # mixed up order is intentional
    C = ['sample_data/sample_a.bam']
    D = ['sample_data/sample_b.bam']
    B = ['sample_data/sample_b.bam']
",
            test_path
        );
        use std::path::Path;
        let output_dir = Path::new(test_path);
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        snp_diff_from_toml(&toml).unwrap();
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
