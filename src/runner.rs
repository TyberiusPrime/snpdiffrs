#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unused_imports)]
use crate::chunked_genome;
use crate::chunked_genome::Chunk;
use crate::consts::*;
use crate::coverage::{Coverage, ResultRow};
use itertools::{iproduct, Itertools};
use rust_htslib::bam;
use serde::Deserialize;
use slog::info;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::{RwLock, RwLockWriteGuard};
use toml;

use crossbeam::crossbeam_channel::unbounded;
use crossbeam::thread;

fn default_quality_threshold() -> u8 {
    15u8
}

fn default_blocksize() -> u32 {
    50_000_000
}

#[derive(Deserialize, Debug)]
pub struct RunConfig {
    output_dir: String,
    chromosomes: Option<Vec<String>>,
    samples: HashMap<String, Vec<String>>,
    #[serde(default = "default_blocksize")]
    block_size: u32,
    ncores: Option<u32>,
    #[serde(default = "default_quality_threshold")]
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: Option<f32>,
}

impl RunConfig {
    fn input_filenames(&self) -> Vec<Vec<PathBuf>> {
        self.samples
            .iter()
            .map(|(name, bam_files)| bam_files.iter().map(|x| x.into()).collect())
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
        .into_iter()
        .next()
        .expect("No samples")
        .iter()
        .next()
        .unwrap();
    let first_bam = bam::IndexedReader::from_path(first_bam).unwrap();

    let pairs: Vec<_> = config.samples.keys().sorted().combinations(2).collect();
    {
       let chunks = chunked_genome::ChunkedGenome::new(first_bam, &config.chromosomes);
        let mut runner = NtoNRunner::new(config, chunks.iter(50_000_000).collect());
        runner.run();
    }
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

enum BlockStatus {
    Unloaded,
    Loaded(usize, usize, Coverage), //for which chunk, (input is by position in block_states), remaining usages, coverage
    ReloadPlease(usize),
    Done,
}
/*

struct Runner {
    config: RunConfig,
    input_filenames: Vec<Vec<PathBuf>>,
    output_handles: Vec<std::io::BufWriter<std::fs::File>>, //one output handle per pair
    pairs: Vec<(usize, usize)>,
    chunks: Vec<chunked_genome::Chunk>,
    block_unload_threshold: usize,
    n_threads: usize,


    block_states: Vec<BlockStatus>,
    pairs_remaining: Vec<(usize, (usize, usize))>, // chunk, pair
}

impl Runner {
    fn new(
        config: RunConfig,
        output_handles: Vec<std::io::BufWriter<std::fs::File>>, //one output handle per pair
        pairs: Vec<(usize, usize)>,
        chunks: Vec<chunked_genome::Chunk>,
        n_threads: usize,
    ) -> Runner {
        let chunk_size = (chunks[0].stop - chunks[0].start) as usize;
        let input_filenames = config.input_filenames();
        let block_unload_threshold: usize = if pairs.len() == input_filenames.len() - 1 {
            input_filenames.len() - 1
        } else if pairs.len() == (input_filenames.len().pow(2) as f64 / 2.).ceil() as usize {
            (input_filenames.len().pow(2) as f64 / 2.).ceil() as usize
        } else {
            panic!("Could not calculate block_unload_threshold - bug");
        };
        let block_count = input_filenames.len();
        let mut block_states = Vec::new();
        for ii in 0..block_count {
            block_states.push(BlockStatus::Unloaded);
        }


        let mut pairs_remaining: Vec<(usize, (usize, usize))> = Vec::new();
        for ii in 0..chunks.len() {
            for p in pairs.iter() {
                pairs_remaining.push((ii, *p));
            }
        }

        Runner{
            config,
            input_filenames,
            output_handles,
            pairs,
            chunks,
            block_unload_threshold,
            n_threads,
            block_states,
            pairs_remaining,
        }
    }

    fn run(mut self){
        use threadpool::ThreadPool;
        let pool = ThreadPool::new(self.n_threads);
        let rw_block_states = RwLock::new(self.block_states);
        let rw_pairs_remaining = RwLock::new(self.pairs_remaining);
        let input_filenames = self.input_filenames;
        let quality_threshold = self.config.quality_threshold;
        let max_chunk = self.chunks.len() - 1;
        let chunks = self.chunks;
        let filter_homo_polymer_threshold = self.config.filter_homo_polymer_threshold.clone();
        let pairs_per_block = self.block_unload_threshold;
        use slog::{Drain, o, info};

        let decorator = slog_term::TermDecorator::new().build();
        let drain = slog_term::CompactFormat::new(decorator).build().fuse();
        let drain = slog_async::Async::new(drain).build().fuse();

        let log = slog::Logger::root(drain, o!());

        info!(log, "Logging ready!");
        pool.execute(move || {
            loop {
                let input_filenames: Vec<Vec<&Path>> = input_filenames.iter().map(|x| x.iter().map(|y| y.as_path()).collect()).collect();
                let chunks = &chunks;
                let log = log.clone();
                info!(log, "Enter loop");
                let do_load = {Self::unload_blocks(&mut *rw_block_states.write().unwrap(), max_chunk, &log)};
                let mut found = false;
                info!(log, "do load: {:?}", &do_load);
                if let Some((block_no, chunk_no)) = do_load {
                    let chunk = &chunks[chunk_no];
                    let cov = Coverage::from_bams(&input_filenames[block_no], chunk.tid, chunk.start, chunk.stop,
                                                  quality_threshold, &filter_homo_polymer_threshold);
                    rw_block_states.write().unwrap()[block_no] = BlockStatus::Loaded(chunk_no,
                                                                                     pairs_per_block, cov);
                    info!(log, "loaded block_no {}, chunk_no, {}", &block_no, &chunk_no);
                    found = true;
                }
                for (ii, (chunk, pair)) in rw_pairs_remaining.read().unwrap().iter().enumerate() {
                    info!(log, "Examining chunk {}, pair {:?} ", &chunk ,&pair);
                    if let BlockStatus::Loaded(a_chunk_no, _remaining, covA) = &rw_block_states.read().unwrap()[pair.0] {
                    if let BlockStatus::Loaded(b_chunk_no, _remaining, covb) = &rw_block_states.read().unwrap()[pair.1] {
                        info!(log, "Both blocks loaded");
                        if *a_chunk_no == *chunk && *b_chunk_no == *chunk {
                            info!(log, "Calculating chunk: {}, pair: {:?}", &chunk ,&pair);
                            &rw_pairs_remaining.write().unwrap().remove(ii);
                            info!(log, "poped from todo");
                            //insert execution and calculation here.
                            {
                                let bs = &mut rw_block_states.write().unwrap();
                                if let BlockStatus::Loaded(_, remaining, _) = &mut bs[pair.0] {
                                    *remaining -= 1;
                                }
                                if let BlockStatus::Loaded(_, remaining, _) = &mut bs[pair.1] {
                                    *remaining -= 1;
                                }
                            }
                            found = true;
                            break;
                        } else
                        {
                            info!(log, "Wrong chunks");
                        }
                    }
                    }
                }
                if !found {
                    info!(log, "Had to sleep");
                    std::thread::sleep(std::time::Duration::new(1,0)); //todo: wake on semaphore taht a block was loaded?
                }


            }

        });
        pool.join();
        }

    fn unload_blocks(block_states: &mut Vec<BlockStatus>, max_chunk: usize, log: &slog::Logger) -> Option<(usize, usize)> {
        let mut do_load = None;
        for (ii, status) in block_states.iter_mut().enumerate() {
            match status {
                BlockStatus::Loaded(chunk, remaining, _)=> {
                    let next_chunk = *chunk + 1;
                    if *remaining == 0 {
                        if next_chunk < max_chunk {
                            info!(&log, "unloading {}, chunk: {}", &ii, &chunk);
                            *status = BlockStatus::ReloadPlease(next_chunk);
                            do_load = Some((ii, next_chunk));
                            break;
                        } else {
                            *status = BlockStatus::Done
                        }
                    };
                    },
                BlockStatus::Unloaded => {
                    *status = BlockStatus::ReloadPlease(0);
                    do_load = Some((ii, 0));
                    break;
                }
                _ => {}
            }

        }
        do_load
    }
    */

fn get_logger() -> slog::Logger {
    use slog::{info, o, Drain};

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    slog::Logger::root(drain, o!())
}

enum ToDo {
    LoadCoverage(Arc<Vec<PathBuf>>, usize, Chunk, usize),
    CalcSnps(usize, [usize; 2], Arc<Coverage>, Arc<Coverage>, Chunk),
    OutputResult(Vec<ResultRow>, Chunk, Arc<Mutex<BufWriter<File>>>),
    Quit,
}

impl std::fmt::Debug for ToDo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ToDo::LoadCoverage(_, _block_no, _, _chunk_no) => {
                f.write_fmt(format_args!("ToDo::LoadCoverage"))
            }
            ToDo::CalcSnps(_, _, _, _, _) => f.write_fmt(format_args!("ToDo:Calc_Snps")),
            ToDo::OutputResult(_,_, _) => f.write_fmt(format_args!("ToDo:OutputResult")),
            ToDo::Quit => f.write_fmt(format_args!("ToDo:Quit")),
        }
    }
}

enum JobResult {
    LoadedCoverage(usize, usize, Coverage, Chunk), //chunk, input_id,
    CalculatedSnps(usize, [usize; 2], Vec<ResultRow>, Chunk),
    QuitDone,
}

impl std::fmt::Debug for JobResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            JobResult::LoadedCoverage(chunk_id, input_id, _, _chunk_no) => {
                f.write_fmt(format_args!(
                    "JobResult::LoadedCoverage(chunk: {}, input_id: {})",
                    chunk_id, input_id
                ))
            }
            JobResult::CalculatedSnps(chunk_id, pair, result_rows, _chunk) => f.write_fmt(format_args!(
                "JobResult::CalculatedSnp(c: {}, pair: {:?}, len(results)=={}",
                chunk_id,
                pair,
                result_rows.len()
            )),
            JobResult::QuitDone => f.write_fmt(format_args!("JobResult::QuitDone")),
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
    input_id: usize, 
    coverage: Arc<Coverage>,
    used: usize,
}

impl NtoNRunner {
    fn new(config: RunConfig, chunks: Vec<Chunk>) -> NtoNRunner {
        NtoNRunner { config, chunks }
    }

    fn run(self) {
        let log = get_logger();

        let mut blocks: Vec<Block> = Vec::new(); //chunk_id, input_id, coverage, usage_count
        let input_files: Vec<Vec<PathBuf>> = self
            .config
            .input_filenames()
            .iter()
            .map(|x| x.iter().map(|y| PathBuf::from(y)).collect())
            .collect();
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_files.into_iter().map(|x| Arc::new(x)).collect();
        let block_iterator: Vec<_> =
            iproduct!(self.chunks.iter().enumerate(), 0..input_filenames.len())
                .enumerate()
                .collect();
        let mut block_iterator: Vec<_> = block_iterator.into_iter().rev().collect();
        let quality_threshold = self.config.quality_threshold;
        let filter_homo_polymer_threshold = self.config.filter_homo_polymer_threshold.clone();
        let min_score = self.config.min_score.unwrap_or(50.0);

        let samples: Vec<&String> = self.config.samples.keys().collect();

        let output_dir = Path::new(&self.config.output_dir);
        let pairs: Vec<[usize; 2]> = (0..input_filenames.len())
            .combinations(2)
            .map(|x| [x[0], x[1]])
            .collect();
        let output_filenames: Vec<_> = pairs
            .iter()
            .map(|x| output_dir.join(format!("{}_vs_{}.tsv.tmp", samples[0], samples[1])))
            .collect();
        let outputs: Vec<_> = output_filenames
            .iter()
            .map(|filename| std::fs::File::create(filename).unwrap())
            .map(|handle| std::io::BufWriter::new(handle))
            .map(|buf| Arc::new(Mutex::new(buf)))
            .collect();

        let mut outputs: HashMap<[usize; 2], _> = pairs.iter().map(|x| *x).zip(outputs).collect();

        let ncores = self.config.ncores.unwrap_or(num_cpus::get() as u32);
        let max_concurrent_blocks = input_filenames.len() + ncores as usize; //we need this many, right?
        info!(log, "ncores: {}", ncores);

        let (todo_sender, todo_receiver) = crossbeam::crossbeam_channel::unbounded();
        for ii in 0..max_concurrent_blocks {
            if block_iterator.is_empty() {
                if ii == 0 {
                    panic!("Nothing to do?");
                }
                else  {
                    break;
                }
            }
            let (_, ((chunk_no, chunk), block_no)) = block_iterator.pop().unwrap();
            todo_sender
                .send(ToDo::LoadCoverage(
                    input_filenames[block_no].clone(),
                    block_no,
                    chunk.clone(),
                    chunk_no,
                ))
                .expect("Could not send initial loads");
        }

        let (result_sender, result_receiver) = crossbeam::crossbeam_channel::unbounded();
        let mut block_iterator = Arc::new(Mutex::new(block_iterator));
        let mut quit_counter = 0;
        thread::scope(|s| {
            let mut block_iterator = block_iterator.clone();
            for nn in 0..ncores {
                let thread_recv = todo_receiver.clone();
                let result_sender = result_sender.clone();
                let log = log.clone();
                s.spawn(move |s2| {
                    //so this is what happens in the worker threads.
                    for el in thread_recv {
                        match el {
                            ToDo::LoadCoverage(input_filenames, input_no, chunk, chunk_no) => {
                                info!(log, "exc: Load_coverage: c:{} i:{}", chunk_no, input_no);
                                result_sender
                                    .send(JobResult::LoadedCoverage(
                                        chunk_no,
                                        input_no,
                                        Coverage::from_bams(&input_filenames.as_ref().iter().map(|x| x.as_path()).collect(),
                                                            chunk.tid,
                                                            chunk.start,
                                                            chunk.stop,
                                                            quality_threshold,
                                                            &filter_homo_polymer_threshold,
                                                            ) ,
                                        chunk,
                                    ))
                                    .expect("Could not send LoadedCoverage reply");
                            }
                            ToDo::CalcSnps(chunk_id, pair, cov, other_cov, chunk) => {
                                let result: Vec<ResultRow> = other_cov.score_differences(cov.as_ref(), min_score, chunk.start);
                                info!(log, "exc: calc snps c: {}, pair: {:?}, result_len: {}", chunk_id, pair, result.len());
                                result_sender.send(
                                    JobResult::CalculatedSnps(chunk_id, pair, result, chunk)).expect("Could not send CalculatedSnps");
                            }
                            ToDo::OutputResult(result_rows, chunk, output) => {
                                info!(log, "exc: output results for chr: {}, start: {}, count: {}, first entry: {}", 
                                      chunk.chr, chunk.start, result_rows.len(),
                                      if result_rows.is_empty() {0} else {result_rows[0].relative_pos});
                                write_results(
                                    &mut output.lock().unwrap(),
                                    &chunk.chr, 
                                    0, result_rows);
                            }
                            ToDo::Quit => {
                                info!(log, "exc: quit: {}", nn);
                                result_sender.send(JobResult::QuitDone).expect("Could not send QuitDone");
                                return;
                            }
                        }
                    }
                });
            }

            for res in result_receiver {
                info!(log, "Received a result: {:?}", res);
                match res {
                    JobResult::LoadedCoverage(chunk_id, input_id, cov, chunk) => {
                        //Todo: optimize by not putting empty blocks into the pool
                        let cov = Arc::new(cov);
                        let mut used = 0;
                        if blocks.len() > 0 {
                            for block in
                                blocks.iter_mut()
                            {
                                info!(log, "checking block other_chunk_id {}: other_input_id, {}", block.chunk_id, block.input_id);
                                if block.chunk_id == chunk_id{
                                    let (a,b, cov_a, cov_b) =
                                        if block.input_id < input_id {
                                            (block.input_id, input_id, cov.clone(), block.coverage.clone())
                                        }
                                        else  {
                                            (input_id, block.input_id, block.coverage.clone(), cov.clone())
                                        };
                                    info!(
                                        log,
                                        "Sending pair {} {} for chunk {}",
                                        a,
                                        b,
                                        chunk_id
                                        );
                                    todo_sender
                                        .send(ToDo::CalcSnps(
                                                chunk_id,
                                                [a, b],
                                                cov_a,
                                                cov_b,
                                                chunk.clone(),
                                                ))
                                        .expect("Could not send CalcSnps");
                                    block.used += 1;
                                    used += 1;
                                }
                            }
                        }
                        blocks.push(Block{chunk_id, input_id, coverage: cov, used});
                        },
                    JobResult::CalculatedSnps(chunk_id, pair, result_rows, chunk) =>{
                        info!(log, "Received CalculatedSnps for c:{}, pair: {:?}", chunk_id, pair);
                        todo_sender.send(ToDo::OutputResult(result_rows, chunk,
                                                            outputs.get(&pair).expect("Could not find output for pair").clone())
                                ).expect("Could not send OutputResult");
                        let mut to_delete = Vec::new();
                        let mut found = 0;
                        for (ii, block) in blocks.iter_mut().enumerate() {
                            if (block.chunk_id == chunk_id) && ((block.input_id == pair[0]) || (block.input_id == pair[1])) {
                                info!(log, "Accepted block: {:?}", block);
                                block.used -= 1;
                                found += 1;
                                if block.used == 0{
                                    info!(log, "deleting block c: {}, i: {}; ii: {}", block.chunk_id, block.input_id, ii);
                                    to_delete.push(ii);
                                }
                            }
                        }
                        if found != 2 {
                            panic!("you made a mistake and we got back a block that no longer exists: c: {}, pair: {:?}, found: {}", chunk_id, pair, found);
                        }
                        for ii in to_delete.iter().rev() {
                            info!(log, "deleting block at {}", ii);
                            blocks.remove(*ii);
                        }
                        if block_iterator.lock().unwrap().len() == 0 { //must do this first, not after sending out some more work...
                            info!(log, "block_iterator empty - preparing to leave");
                            for nn in 0..ncores {
                                todo_sender.send(ToDo::Quit).unwrap();
                            }
                        }

                        if blocks.len() < max_concurrent_blocks {
                            for cc in blocks.len()..max_concurrent_blocks {
                                let next = block_iterator.lock().expect("Could not lock block_iterator").pop(); 
                                if let Some((_, ((chunk_no, chunk), block_no))) = next {
                                        todo_sender
                                            .send(ToDo::LoadCoverage(
                                                input_filenames[block_no].clone(),
                                                block_no,
                                                chunk.clone(),
                                                chunk_no,
                                            ))
                                            .unwrap();
                                }
                            }
                        }
                    },
                    JobResult::QuitDone => {
                        quit_counter += 1;
                        if quit_counter == ncores {
                            break;
                        }

                    }
                }
            }
            info!(log, "Left receive results loop");
        })
        .unwrap();

        output_filenames
            .iter()
            .map(|tmp_filename| {
                std::fs::rename(
                    tmp_filename,
                    tmp_filename.with_file_name(tmp_filename.file_stem().unwrap()),
                )
                .unwrap()
            })
            .count();
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
}
