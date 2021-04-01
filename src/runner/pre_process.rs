use super::{all_files_exists, ensure_output_dir, get_logger, haplotype_const_to_str, RunConfig};
use crate::chunked_genome::{Chunk, ChunkedGenome};
use crate::coverage::{Coverage, EncCoverage, ResultRow};

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
    all_files_exists(&config.samples, "samples");
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

pub fn name_folder(
    sample_name: &str,
    input_filenames: &Vec<PathBuf>,
    block_size: usize,
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
) -> String {
    use sha2::{Digest, Sha256};
    use std::os::unix::ffi::OsStrExt;
    let mut hasher = Sha256::new();
    for p in input_filenames {
        hasher.update(p.as_os_str().as_bytes())
    }
    let hash = hasher.finalize();
    format!(
        "{}_bs={}_q={}_homo_polymer={}_{}",
        sample_name,
        block_size,
        quality_threshold,
        match filter_homo_polymer_threshold {
            Some(x) => format!("{}", x),
            None => "None".to_string(),
        },
        format!("{:x}", hash)
    )
    .to_string()
}

pub fn name_block(chr: &str, start: u32) -> String {
    format!("{}_{}.snpdiffrs_block", chr, start)
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
            .map(|(sample_name, input_files)| {
                let od = output_dir.join(name_folder(
                    sample_name,
                    input_files,
                    config.block_size,
                    config.quality_threshold,
                    config.filter_homo_polymer_threshold,
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
                            od.join(name_block(&chunk.chr, chunk.start)),
                        )
                    })
            })
            .flatten()
            .collect();

        let _: Vec<()> = todo
            .par_iter()
            .map(|(input_filenames, chunk, output_filename)| {
                if !output_filename.exists() {
                    //note that no coverage regions get to be cov of 0 bytes
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
                    let cov = EncCoverage::new(&cov);
                    cov.to_file(output_filename)
                        .expect("Could not write to file");
                }
            })
            .collect();
    }
}
