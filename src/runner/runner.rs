use std::sync::Arc;
use std::sync::Mutex;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};


#[allow(unused_imports)]
use slog::{debug, info};
use itertools::{iproduct, Itertools};
use crossbeam::thread;

use crate::coverage::{Coverage, EncCoverage, ResultRow};
use crate::chunked_genome::{Chunk, ChunkedGenome};
use super::{all_files_exists, ensure_output_dir, get_logger, haplotype_const_to_str, RunConfig};

pub trait RunStrategy: Send + Sync{
    fn expected_number_of_comparisons(&self) -> usize;
    fn expected_number_of_uses_per_block(&self, sample_name: &str) -> usize;
    fn include_comparison(&self, sample_a_name: &str, sample_b_name: &str) -> bool;
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
    coverage_a: Arc<EncCoverage>,
    coverage_b: Arc<EncCoverage>,
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
    coverage: EncCoverage,
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
pub struct Runner {
    //config: RunConfig,
    //chunks: Vec<Chunk>,
    strategy: Arc<Box<dyn RunStrategy>>,
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
    coverage: Arc<EncCoverage>,
    remaining_uses: usize,
}

impl Runner 
{
    pub fn new(config: RunConfig, samples: Vec<String>, input_filenames: Vec<Vec<PathBuf>>, run_strategy: Box<dyn RunStrategy>) -> Runner {
        let log = get_logger();
        //input
        if samples.len() < 2 {
            panic!("Not enough samples supplied");
        }
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_filenames.into_iter().map(Arc::new).collect();

        //output
        //this is only right for the n-to-n case...
        let output_dir = Path::new(&config.output_dir);
        let mut outputs: HashMap<[usize; 2], _> = HashMap::new();
        let mut output_filenames = Vec::new();
        for pair in (0..input_filenames.len()).combinations(2) { //the samples are sorted by sample name.
            let sample_a = &samples[pair[0]];
            let sample_b = &samples[pair[1]];
            if run_strategy.include_comparison(sample_a, sample_b) {
                let filename = output_dir.join(format!("{}_vs_{}.tsv.tmp", sample_a, sample_b));
                let fh = std::fs::File::create(&filename).expect("Could not create output file");
                let fh = Arc::new(Mutex::new(std::io::BufWriter::new(fh)));
                outputs.insert([pair[0], pair[1]], fh);
                debug!(log, "output for {}{} is {:?}", pair[0], pair[1], &filename);
                output_filenames.push(filename);
            }
        }
        debug!(log, "outputs: {}, expected comparisons: {}", outputs.len(), run_strategy.expected_number_of_comparisons());
        assert!(outputs.len() == run_strategy.expected_number_of_comparisons());



        let ncores = config.ncores.unwrap_or(num_cpus::get() as u32);
        let max_concurrent_blocks =
            (input_filenames.len() + ncores as usize).max(ncores as usize * 3); //we need at least one per input_filename to create all pairs. And some bonus is helpful to avoid delays.
        debug!(log, "max concurrent blocks {}", max_concurrent_blocks);

        let chunks: Vec<Chunk> = ChunkedGenome::new(config.first_bam(), &config.chromosomes).iter(config.block_size).collect();

        Runner {
            strategy: Arc::new(run_strategy),
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

    pub fn run(self) {
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

        let expected_number_of_outputs = self.strategy.expected_number_of_comparisons() * chunk_count;

        /*
        let mut load_progress = MappingBar::with_range(
            0,
            block_iterator.len()
                + block_iterator.len() * expected_number_of_uses_per_block * chunk_count,
        )
        .timed();
        load_progress.set_len(40);
        load_progress.set(0usize);
        print!("{}", load_progress);
        */
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
            "Expected number of output blocks: {}", expected_number_of_outputs
        );
        let mut blocks: Vec<Block> = Vec::new(); //the currently loaded blocks. I suppose those should no texceed max_concurrent_blocks
        let remaining_outputs = std::sync::atomic::AtomicUsize::new(expected_number_of_outputs);

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
                            let cov: Option<EncCoverage> = match _self.preprocessed_dir.as_ref() {
                                Some(pd) => {
                                    let full_pd = pd.join(super::pre_process::name_folder(
                                        &_self.samples[payload.sample_id],
                                        &_self.input_filenames[payload.sample_id],
                                        _self.block_size,
                                        _self.quality_threshold,
                                        _self.filter_homo_polymer_threshold));
                                    let filename = full_pd.join(super::pre_process::name_block(&payload.chunk.chr, payload.chunk.start));
                                    if !filename.exists() {
                                        //info!(log, "Failed to load preprocessed from: {:?}", filename);
                                        None
                                    } else {
                                        EncCoverage::from_file(&filename)
                                    }
                                },
                                None => None
                            };
                            let cov: EncCoverage = match cov {
                                Some(cov) => cov,
                                None => {
                                        let e = EncCoverage::new(&Coverage::from_bams(&_self.input_filenames[payload.sample_id].as_ref().iter().map(|x| x.as_path()).collect::<Vec<_>>(),
                                        payload.chunk.tid,
                                        payload.chunk.start,
                                        payload.chunk.stop,
                                        _self.quality_threshold,
                                        &_self.filter_homo_polymer_threshold,
                                        ));
                                        if let Some(pd) = _self.preprocessed_dir.as_ref() {
                                            let full_pd = pd.join(super::pre_process::name_folder(
                                                &_self.samples[payload.sample_id],
                                                &_self.input_filenames[payload.sample_id],
                                                _self.block_size,
                                                _self.quality_threshold,
                                                _self.filter_homo_polymer_threshold));
                                            if !full_pd.exists() {
                                                std::fs::create_dir_all(&full_pd).unwrap();
                                            }
                                            let filename = full_pd.join(super::pre_process::name_block(&payload.chunk.chr, payload.chunk.start));
                                            info!(log, "Rebuilding preprocessed block: {:?}", filename);
                                            //with writing: 45.19
                                            //without: 26s.
                                            //after writing: 3s
                                            e.to_file(&filename).expect("Could not write preprocessed file");
                                        }
                                        e
                                }
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
                            debug!(log, "Calced {}, {} for {} {}" ,payload.chunk.chr, payload.chunk.start, payload.sample_pair[0], payload.sample_pair[1]);
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
                        debug!(log, "loaded {} {} {} {}", payload.chunk.chr, payload.chunk.start, payload.coverage.len(), payload.coverage.sum());
                        /*load_progress.add(1usize);
                        if load_progress.has_progressed_significantly() {
                            print!("\r{}", load_progress);
                        }
                        */

                        let cov = Arc::new(payload.coverage);
                        let remaining_uses = _self.strategy.expected_number_of_uses_per_block(&_self.samples[payload.sample_id]);
                        if !blocks.is_empty() {
                            let block_count = blocks.len(); //debug
                            for block in
                                blocks.iter_mut()
                            {
                                if block.chunk_id == payload.chunk_id{
                                    let (a,b, cov_a, cov_b) =
                                        //do I need the string comparison here?
                                        if _self.samples[block.sample_id] < _self.samples[payload.sample_id] {
                                            (block.sample_id, payload.sample_id, cov.clone(), block.coverage.clone())
                                        }
                                        else  {
                                            (payload.sample_id, block.sample_id, block.coverage.clone(), cov.clone())
                                        };
                                    if _self.strategy.include_comparison(&_self.samples[a], &_self.samples[b]) {
                                        debug!(log, "Sending pair {} {} for chunk {} {}. block count: {}", _self.samples[a], _self.samples[b], payload.chunk.chr, payload.chunk.start, block_count);
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
                                    } else {
                                        debug!(log, "Skipped because of strategy: pair {} {} for chunk {} {}. block count: {}", _self.samples[a], _self.samples[b], payload.chunk.chr, payload.chunk.start, block_count);
                                    }
                                }
                            }
                        }
                        debug!(log,"Pushing block c: {}, sample: {}, remaining: {}", payload.chunk_id, _self.samples[payload.sample_id], remaining_uses);
                        blocks.push(Block{chunk_id: payload.chunk_id, sample_id: payload.sample_id, coverage: cov, remaining_uses});
                    },

                    MsgDone::CalcSnps(payload) => {
                            debug!(log, "Received CalcSnps result for c:{}, pair: {:?}", payload.chunk_id, payload.sample_pair);
                            /*
                            load_progress.add(1usize);
                            if load_progress.has_progressed_significantly() {
                                print!("\r{}", load_progress);
                            }
                            */
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
                            let deleted_count = to_delete.len();
                            for ii in to_delete.iter().rev() {
                                debug!(log, "deleting block at {} - arc count strong: {} weak; {}", ii, 
                                    Arc::strong_count(&blocks[*ii].coverage),
                                    Arc::weak_count(&blocks[*ii].coverage),
                                );
                                blocks.remove(*ii);
                            }
                            if deleted_count > 0 {
                                let mut bi = block_iterator.lock().unwrap(); //let's lock this for the whole block
                                debug!(log, "Pushing {} new blocks", deleted_count);
                                for _ in 0..deleted_count {
                                    Self::_push_next_block(&mut bi, &todo_sender);
                                }
                            }
                        },
                    MsgDone::OutputDone => {
                            let rr = remaining_outputs.fetch_sub(1, std::sync::atomic::Ordering::SeqCst);
                            debug!(log, "MsgDone:OutputDone - remaining outputs {}", rr);
                            if rr == 1 { //previous value!
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
