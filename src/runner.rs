use crate::chunked_genome;
use crate::consts::*;
use crate::coverage::{Coverage, ResultRow};
use itertools::Itertools;
use rust_htslib::bam;
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;
use toml;

fn default_quality_threshold() -> u8 {
    15u8
}

#[derive(Deserialize, Debug)]
pub struct RunConfig {
    output_dir: String,
    chromosomes: Option<Vec<String>>,
    samples: HashMap<String, Vec<String>>,
    #[serde(default = "default_quality_threshold")]
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: Option<f32>,
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
    let output_filenames: Vec<_> = pairs
        .iter()
        .map(|x| output_dir.join(format!("{}_vs_{}.tsv.tmp", x[0], x[1])))
        .collect();
    {
        let mut output_handles: Vec<_> = output_filenames
            .iter()
            .map(|filename| std::fs::File::create(filename).unwrap())
            .map(|handle| std::io::BufWriter::new(handle))
            .collect();

        let quality_threshold = config.quality_threshold;
        let filter_homo_polymer_threshold = config.filter_homo_polymer_threshold.clone();
        let min_score = config.min_score.unwrap_or(50.0);

        let chunks = chunked_genome::ChunkedGenome::new(first_bam, &config.chromosomes);
        for chunk in chunks.iter(50_000_000) {
            let cov = get_coverages(
                &config,
                &chunk,
                quality_threshold,
                &filter_homo_polymer_threshold,
            );
            calculate_differences(&cov, &pairs, &mut output_handles, &chunk, min_score);
        }
    }
    //files are now closed. rename them all
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

    Ok(())
}

fn get_coverages(
    config: &RunConfig,
    chunk: &chunked_genome::Chunk,
    quality_threshold: u8,
    filter_homo_polymer_threshold: &Option<u8>,
) -> HashMap<String, Coverage> {
    config
        .samples
        .iter()
        .map(|(name, bam_files)| {
            (
                name.clone(),
                Coverage::from_bams(
                    bam_files,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    quality_threshold,
                    filter_homo_polymer_threshold,
                ),
            )
        })
        .collect()
}

fn calculate_differences(
    coverages: &HashMap<String, Coverage>,
    pairs: &Vec<Vec<&String>>,
    output_handles: &mut Vec<std::io::BufWriter<std::fs::File>>,
    chunk: &chunked_genome::Chunk,
    min_score: f32,
) {
    for (ii, pair) in pairs.iter().enumerate() {
        let delta = coverages[pair[0]].score_differences(&coverages[pair[1]], min_score);
        if !delta.is_empty() {
            write_results(&mut output_handles[ii], &chunk.chr, chunk.start, delta);
        }
    }
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
