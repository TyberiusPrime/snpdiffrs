use super::{all_files_exists, ensure_output_dir, get_logger, haplotype_const_to_str, RunConfig};
use crate::chunked_genome::{Chunk, ChunkedGenome};
use crate::coverage::{Coverage, ResultRow};

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use lzzzz::{lz4, lz4_hc, lz4f};
use rayon::prelude::*;

/* the preprocessed format is simple.
 * You create one folder per sample/block_size, and one file per chromosome/start_pos
 * Note that the block_size is constant within one run.
 *
*/

pub fn run_preprocess(config: RunConfig) -> Result<(), ()> {
    all_files_exists(&config.samples);
    ensure_output_dir(&Path::new(&config.output_dir));

    let chunks = ChunkedGenome::new(config.first_bam(), &config.chromosomes);
    let block_size = config.block_size;
    let runner = PreProcessRunner::new(config, chunks);
    runner.run();
    Ok(())
}

struct PreProcessRunner {
    input_filenames: Vec<Arc<Vec<PathBuf>>>,
    output_directories: Vec<PathBuf>,
    chunks: Vec<Chunk>,
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
}

impl PreProcessRunner {
    fn new(config: RunConfig, chunked_genome: ChunkedGenome) -> PreProcessRunner {
        //input
        let samples_and_input_filenames = config.samples_and_input_filenames();
        let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) =
            samples_and_input_filenames.clone().into_iter().unzip();
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_filenames.into_iter().map(Arc::new).collect();
        //output
        let output_dir = Path::new(&config.output_dir);
        let output_directories = samples_and_input_filenames
            .iter()
            .map(|(sample_name, _input_files)| {
                let od = output_dir.join(format!(
                    "{}_bs={}_q={}_homo_polymer={}.snpdiff_block",
                    sample_name,
                    config.block_size,
                    config.quality_threshold,
                    match config.filter_homo_polymer_threshold {
                        Some(x) => format!("{}", x),
                        None => "None".to_string(),
                    }
                ));
                if !od.exists() {
                    std::fs::create_dir_all(&od).unwrap();
                }
                od
            })
            .collect();

        let chunks_vec: Vec<Chunk> = chunked_genome.iter(config.block_size).collect();
        PreProcessRunner {
            input_filenames,
            output_directories,
            chunks: chunks_vec,
            quality_threshold: config.quality_threshold,
            filter_homo_polymer_threshold: config.filter_homo_polymer_threshold,
        }
    }

    fn run(self) {
        let todo: Vec<(Arc<Vec<PathBuf>>, &Chunk, PathBuf)> = self
            .chunks
            .iter()
            .map(|chunk| {
                self.input_filenames
                    .iter()
                    .zip(self.output_directories.iter())
                    .map(move |(inputs, od)| {
                        (
                            inputs.clone(),
                            chunk,
                            od.join(
                                format!("{}_{}.snpdiffrs_block", chunk.chr, chunk.start).to_owned(),
                            ),
                        )
                    })
            })
            .flatten()
            .collect();

        let _: Vec<()> = todo
            .par_iter()
            .map(|(input_filenames, chunk, output_filename)| {
                if !output_filename.exists() {
                    //todo: remove low coverage regions for smaller files
                    let cov = Coverage::from_bams(
                        &input_filenames
                            .as_ref()
                            .iter()
                            .map(|x| x.as_path())
                            .collect::<Vec<_>>(),
                        chunk.tid,
                        chunk.start,
                        chunk.stop,
                        self.quality_threshold,
                        &self.filter_homo_polymer_threshold,
                    );

                    //lz4_hc::compress_to_vec(data, &mut comp, lz4_hc::CLEVEL_DEFAULT)?;
                    let mut fh = File::create(output_filename.clone())
                        .expect(&format!("Unable to create file {:?}", output_filename));
                    use safe_transmute::{transmute_to_bytes, SingleManyGuard};

                    let raw = cov.into_raw_vec();
                    let raw_v8: &[u8] = transmute_to_bytes(&raw[..]);
                    //zstd::stream::copy_encode(raw_v8, fh, 9).unwrap();
                    fh.write_all(raw_v8).unwrap();
                }
            })
            .collect();
    }
}
