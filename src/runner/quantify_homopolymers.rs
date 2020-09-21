use super::RunConfig;
use crate::chunked_genome::{Chunk, ChunkedGenome};
use crate::coverage::{count_homopolymers_in_bam, BaseCounts};
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::Mutex;

pub fn run_quantify_homopolymers(config: RunConfig) -> Result<(), ()> {
    for filenames in config
        .samples
        .as_ref()
        .expect("No samples provided")
        .values()
    {
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
    let chunks = ChunkedGenome::new(config.first_bam(), &config.chromosomes);
    let block_size = config.block_size;
    let runner = QuantifyHomoPolymersRunner::new(config, chunks.iter(block_size).collect());
    runner.run();
    Ok(())
}

struct QuantifyHomoPolymersRunner {
    sample_names: Vec<String>,
    input_filenames: Vec<Arc<Vec<PathBuf>>>,
    output_filenames: Vec<PathBuf>,
    outputs: Vec<std::io::BufWriter<std::fs::File>>,
    ncores: u32,
    chunks: Vec<Chunk>,
}

impl QuantifyHomoPolymersRunner {
    fn new(config: RunConfig, chunks: Vec<Chunk>) -> QuantifyHomoPolymersRunner {
        let samples_and_input_filenames = config.samples_and_input_filenames();
        let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) =
            samples_and_input_filenames.into_iter().unzip();
        if samples.len() < 2 {
            panic!("Not enough samples supplied");
        }
        let input_filenames: Vec<Arc<Vec<PathBuf>>> =
            input_filenames.into_iter().map(Arc::new).collect();

        let output_dir = Path::new(&config.output_dir);
        let output_filenames: Vec<_> = samples
            .iter()
            .map(|x| output_dir.join(format!("{}_homopolymer_histogram.tsv.tmp", x)))
            .collect();
        let outputs: Vec<_> = output_filenames
            .iter()
            .map(|filename| std::fs::File::create(filename).unwrap())
            .map(std::io::BufWriter::new)
            //.map(|buf| Arc::new(Mutex::new(buf)))
            .collect();

        let ncores = config.ncores.unwrap_or(num_cpus::get() as u32);

        QuantifyHomoPolymersRunner {
            sample_names: samples,
            input_filenames,
            output_filenames,
            outputs,
            ncores,
            chunks,
        }
    }

    fn run(mut self) {
        let mut counters: Vec<Arc<Mutex<Vec<BaseCounts>>>> = Vec::new();
        for _ in 0..self.input_filenames.len() {
            let v = BaseCounts::new_histogram();
            counters.push(Arc::new(Mutex::new(v)));
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(self.ncores as usize)
            .build_global()
            .unwrap();

        self.chunks.par_iter().for_each(|chunk| {
            for sample in 0..self.input_filenames.len() {
                for filename in self.input_filenames[sample].as_ref().iter() {
                    let here =
                        count_homopolymers_in_bam(filename, chunk.tid, chunk.start, chunk.stop);
                    BaseCounts::update(&mut counters[sample].lock().unwrap(), &here);
                }
            }
        });

        for (sample_name, (out, counter)) in self
            .sample_names
            .iter()
            .zip(self.outputs.iter_mut().zip(counters))
        {
            writeln!(out, "#quantify_homopolymers: {}", sample_name).unwrap();
            writeln!(out, "Homo_polymer\tBase\tCount").unwrap();
            for hp in 0..u8::MAX {
                let hp = hp as usize;
                let row = &counter.lock().unwrap()[hp];
                writeln!(out, "{}\tA\t{}", hp, row.a).unwrap();
                writeln!(out, "{}\tC\t{}", hp, row.c).unwrap();
                writeln!(out, "{}\tG\t{}", hp, row.g).unwrap();
                writeln!(out, "{}\tT\t{}", hp, row.t).unwrap();
            }
        }

        self.output_filenames
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

#[cfg(test)]
mod test {
    use super::super::snp_diff_from_toml;
    use std::path::PathBuf;
    #[test]
    fn test_quantify_homopolymers() {
        let test_dir = "tests/test_quantify_homopolymers";
        let toml = format!(
            "
mode = 'QuantifyHomopolymers'
output_dir = '{}'
block_size = 500000
[samples]
    A_chr4 = ['sample_data/ERR329501_chr4.bam'] # mixed up order is intentional
    B_chr4 = ['sample_data/GSM1553106_chr4.bam'] # mixed up order is intentional
    ",
            test_dir
        );
        use std::path::Path;
        let output_dir = Path::new(test_dir);
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        snp_diff_from_toml(&toml).unwrap();
        let mut tf = PathBuf::new();
        tf.push(output_dir);
        tf.push("A_chr4_homopolymer_histogram.tsv");
        let out: String = std::fs::read_to_string(tf).unwrap();
        assert!(out.contains("\n1\tA\t12\n"));
        assert!(out.contains("\n2\tA\t1468\n"));
        assert!(out.contains("\n2\tC\t1382\n"));
        assert!(out.contains("\n2\tG\t1859\n"));
        assert!(out.contains("\n2\tT\t1149\n"));
    }
}
