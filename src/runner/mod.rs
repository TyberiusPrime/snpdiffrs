use rust_htslib::bam;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

mod n_to_n;
mod pre_process;
mod quantify_homopolymers;

fn default_quality_threshold() -> u8 {
    15u8
}

fn default_blocksize() -> usize {
    50_000_000
}

#[derive(Deserialize, Debug)]
enum RunMode {
    NToN,
    OneToN,
    QuantifyHomopolymers,
    PreProcess,
}

fn default_mode() -> RunMode {
    RunMode::NToN
}

#[derive(Deserialize, Debug)]
pub struct RunConfig {
    #[serde(default = "default_mode")]
    mode: RunMode,
    output_dir: String,
    chromosomes: Option<Vec<String>>,
    samples: Option<HashMap<String, Vec<String>>>,
    references: Option<HashMap<String, Vec<String>>>,
    query: Option<HashMap<String, Vec<String>>>,
    #[serde(default = "default_blocksize")]
    block_size: usize,
    ncores: Option<u32>,
    #[serde(default = "default_quality_threshold")]
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: Option<f32>,
    preprocessed_dir: Option<String>,
}

impl RunConfig {
    fn samples_and_input_filenames(&self) -> Vec<(String, Vec<PathBuf>)> {
        let mut temp: Vec<(&String, &Vec<String>)> = self
            .samples
            .as_ref()
            .expect("No samples provided")
            .iter()
            .collect();
        temp.sort();
        temp.iter()
            .map(|(name, bam_files)| {
                (
                    name.to_string(),
                    bam_files.iter().map(|x| x.into()).collect::<Vec<PathBuf>>(),
                )
            })
            .collect()
    }

    fn first_bam(&self) -> bam::IndexedReader {
        let first_bam = self
            .samples
            .as_ref()
            .expect("No samples provided")
            .values()
            .next()
            .expect("No samples")
            .iter()
            .next()
            .unwrap();
        bam::IndexedReader::from_path(first_bam).unwrap()
    }
}

pub fn snp_diff_from_toml(input: &str) -> Result<(), ()> {
    let config: RunConfig = toml::from_str(input).unwrap();
    match config.mode {
        RunMode::NToN => n_to_n::run_n_to_n(config),
        RunMode::OneToN => panic!("implement me"),
        RunMode::QuantifyHomopolymers => quantify_homopolymers::run_quantify_homopolymers(config),
        RunMode::PreProcess => pre_process::run_preprocess(config),
    }
}

fn haplotype_const_to_str(haplotype: u8) -> &'static str {
    use crate::consts::{AA, AC, AG, AT, CC, CG, CT, GG, GT, NN, TT};
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
pub fn get_logger() -> slog::Logger {
    use slog::{o, Drain};

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).chan_size(20000).build().fuse();

    slog::Logger::root(drain, o!())
}

fn all_files_exists(filenames: &Option<HashMap<String, Vec<String>>>) {
    for sub_filenames in filenames.as_ref().expect("No samples provided").values() {
        for filename in sub_filenames {
            if !Path::new(filename).exists() {
                panic!("File did not exist {}", filename);
            }
        }
    }
}

fn ensure_output_dir(output_dir: &Path) {
    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).unwrap();
    }
}
