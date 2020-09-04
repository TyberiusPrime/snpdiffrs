#![allow(dead_code)]
#![allow(unused_imports)]
use rust_htslib::bam;
use rust_htslib::htslib;
use std::collections::HashMap;
use std::convert::TryInto;
use std::path::Path;
//use rust_htslib::prelude::*;
//#[macro_use] extern crate ndarray;
use ndarray::prelude::*;
use ndarray_stats::QuantileExt;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Read;
use std::io::{Seek, SeekFrom, Write};

#[macro_use]
extern crate lazy_static;
use itertools::Itertools;
use serde::Deserialize;
use toml;

mod chunked_genome;

const BASE_A: usize = 0;
const BASE_C: usize = 1;
const BASE_G: usize = 2;
const BASE_T: usize = 3;
const AA: usize = 0;
const CC: usize = 1;
const GG: usize = 2;
const TT: usize = 3;
const AC: usize = 4;
const AG: usize = 5;
const AT: usize = 6;
const CG: usize = 7;
const CT: usize = 8;
const GT: usize = 9;
const NN: usize = 10;

const READ_ERROR_PROB: f32 = 0.001;
lazy_static! {
    //log functions are not const fn yet
    //these are all named as if the error rate was 1%. (even though it isn't)
    static ref LL_99: f32 = (1.0f32 - READ_ERROR_PROB).ln(); // for AA-haplotypes
    static ref LL_003: f32 = (READ_ERROR_PROB / 3.).ln(); //for AA-haplotypes - misreads
    static ref LL_005: f32 = (READ_ERROR_PROB / 2.).ln(); //for  AC haplotypes - misreads
    static ref LL_495: f32 = ((1. - READ_ERROR_PROB) / 2.).ln(); //AC haplotypes - misreads
    static ref LL_25: f32 = (0.25_f32).ln(); // for nn - not diploid then.
}

fn vector_arg_max(input: &[f32;11]) -> (f32, usize) {
    let mut argmax: usize = 0;
    let mut max = f32::MIN;
    for (ii, value) in input.iter().enumerate() {
        if *value > max {
            argmax = ii;
            max = *value;
        }
    }
    (max, argmax)
}
//
// there is no apperant perfomance difference between u32 + cast and f32 counting
// so we might as well stick to u32 for now.
struct Coverage(Array2<u32>); 

#[derive(Debug)]
struct ResultRow {
    relative_pos: u32,
    count_self_a: u32,
    count_self_c: u32,
    count_self_g: u32,
    count_self_t: u32,
    count_other_a: u32,
    count_other_c: u32,
    count_other_g: u32,
    count_other_t: u32,
    score: f32,
    haplotype_self: u8,
    haplotype_other: u8,
}

impl Coverage {
    fn new(length: usize) -> Coverage {
        Coverage(Array2::zeros((length, 4)))
    }

    fn from_counts(
        count_a: Vec<u32>,
        count_c: Vec<u32>,
        count_g: Vec<u32>,
        count_t: Vec<u32>,
    ) -> Coverage {
        let len = count_a.len();
        if count_c.len() != len || count_g.len() != len || count_t.len() != len {
            panic!("unequal len arrays to from_counts");
        }
        let mut res = Array2::zeros((len, 4));
        for ii in 0..len {
            res[[ii, 0]] = count_a[ii];
            res[[ii, 1]] = count_c[ii];
            res[[ii, 2]] = count_g[ii];
            res[[ii, 3]] = count_t[ii];
        }
        Coverage(res)
    }

    fn from_bam<P: AsRef<Path>>(
        filename: P,
        tid: u32,
        start: u32,
        stop: u32,
        quality_threshold: u8,
        filter_homo_polymer_threshold: &Option<u8>,
    ) -> Coverage {
        Coverage::from_bams(&vec![filename], tid, start, stop, quality_threshold, filter_homo_polymer_threshold)
    }

    fn from_bams<P: AsRef<Path>>(
        bams: &Vec<P>,
        tid: u32,
        start: u32,
        stop: u32,
        quality_threshold: u8,
        filter_homo_polymer_threshold: &Option<u8>,
    ) -> Self {
        let length: u64 = (stop - start).try_into().expect("stop < start");
        let mut result: Coverage = Coverage::new(length as usize);
        let mut any = false;
        for filename in bams {
            any |= result.update_from_bam(
                filename,
                tid,
                start,
                stop,
                quality_threshold,
                filter_homo_polymer_threshold,
            );
        }
        if any {
            result
        } else {
            Coverage::new(0)
        }

    }

    fn update_from_bam<P: AsRef<Path>>(
        &mut self,
        filename: P,
        tid: u32,
        start: u32,
        stop: u32,
        quality_threshold: u8,
        filter_homo_polymer_threshold: &Option<u8>,
    ) -> bool{
        let mut bam = bam::IndexedReader::from_path(filename).expect("Could not read input bam");
        let start = start as i64;
        let stop = stop as i64;
        /*
        let tid: u32 = bam
            .header()
            .tid(chr.as_bytes())
            .expect("Could not find chromosome");
        */
        bam.fetch(tid, start as u64, stop as u64).unwrap();
        let mut read: bam::Record = bam::Record::new();
        let mut any = false;
        while let Ok(true) = bam.read(&mut read) {
            if (read.flags()
                & (htslib::BAM_FUNMAP
                    | htslib::BAM_FSECONDARY
                    | htslib::BAM_FQCFAIL
                    | htslib::BAM_FDUP) as u16)
                > 0
            {
                continue;
            }
            let seq = read.seq();
            if filter_homo_polymer_threshold.is_some()
                && is_homo_polymer(&seq, filter_homo_polymer_threshold.unwrap())
            {
                continue;
            }
            for [read_pos, genome_pos] in read.aligned_pairs().iter() {
                if (*genome_pos < start) || *genome_pos >= stop {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                if read.qual()[*read_pos as usize] < quality_threshold {
                    continue;
                }
                any = true;
                let base: u8 = seq.encoded_base((*read_pos) as usize);
                //base is a 4 bit integer, 0..15 mapping to
                //= 0
                //A 1
                //C 2
                //M 3
                //G 4
                //R 5
                //S 6
                //V 7
                //T 8
                //W 9
                //Y 10
                //H 12
                //K 12
                //D 13
                //B 14
                //N 15
                //no need to reverse, everything is on the same strand...
                let out_base = match base {
                    1 => BASE_A,
                    2 => BASE_C,
                    4 => BASE_G,
                    8 => BASE_T,
                    _ => continue,
                };
                self.0[[(genome_pos - start) as usize, out_base]] += 1;
            }
        }
        any
    }

    fn len(self: &Self) -> usize {
        return self.0.dim().0;
    }

    fn single_log_likelihood(
        count_a: u32,
        count_c: u32,
        count_g: u32,
        count_t: u32,
        which: usize,
    ) -> f32 {
        let a = count_a as f32;
        let c = count_c as f32;
        let g = count_g as f32;
        let t = count_t as f32;
        match which {
            AA => a * *LL_99 + c * *LL_003 + g * *LL_003 + t * *LL_003, //AA
            CC => a * *LL_003 + c * *LL_99 + g * *LL_003 + t * *LL_003, //CC
            GG => a * *LL_003 + c * *LL_003 + g * *LL_99 + t * *LL_003, //GG
            TT => a * *LL_003 + c * *LL_003 + g * *LL_003 + t * *LL_99, //TT
            AC => a * *LL_495 + c * *LL_495 + g * *LL_005 + t * *LL_005, //AC
            AG => a * *LL_495 + c * *LL_005 + g * *LL_495 + t * *LL_005, //AG
            AT => a * *LL_495 + c * *LL_005 + g * *LL_005 + t * *LL_495, //AT
            CG => a * *LL_005 + c * *LL_495 + g * *LL_495 + t * *LL_005, //CG
            CT => a * *LL_005 + c * *LL_495 + g * *LL_005 + t * *LL_495, //CT
            GT => a * *LL_005 + c * *LL_005 + g * *LL_495 + t * *LL_495, //GT
            NN => a * *LL_25 + c * *LL_25 + g * *LL_25 + t * *LL_25,
            _ => panic!("non haplotype passed to single_log_likelihood"),
        }
    }

    fn ll(a: f32, c: f32, g: f32, t:f32) -> [f32;11]{
        let a_ll_003 = a * *LL_003;
        let c_ll_003 = c * *LL_003;
        let g_ll_003 = g * *LL_003;
        let t_ll_003 = t * *LL_003;
        let a_ll_005 = a * *LL_005;
        let c_ll_005 = c * *LL_005;
        let g_ll_005 = g * *LL_005;
        let t_ll_005 = t * *LL_005;
        let a_ll_495 = a * *LL_495;
        let c_ll_495 = c * *LL_495;
        let g_ll_495 = g * *LL_495;
        let t_ll_495 = t * *LL_495;

        [
            a * *LL_99 + c_ll_003 + g_ll_003 + t_ll_003, //AA
            a_ll_003 + c * *LL_99 + g_ll_003 + t_ll_003, //CC
            a_ll_003 + c_ll_003 + g * *LL_99 + t_ll_003, //GG
            a_ll_003 + c_ll_003 + g_ll_003 + t * *LL_99, //TT
            a_ll_495 + c_ll_495 + g_ll_005 + t_ll_005, //AC
            a_ll_495 + c_ll_005 + g_ll_495 + t_ll_005, //AG
            a_ll_495 + c_ll_005 + g_ll_005 + t_ll_495, //AT
            a_ll_005 + c_ll_495 + g_ll_495 + t_ll_005, //CG
            a_ll_005 + c_ll_495 + g_ll_005 + t_ll_495, //CT
            a_ll_005 + c_ll_005 + g_ll_495 + t_ll_495, //GT
            a * *LL_25 + c * *LL_25 + g * *LL_25 + t * *LL_25,
        ]
    }

    fn single_log_likelihood_max_arg_max(
        count_a: u32,
        count_c: u32,
        count_g: u32,
        count_t: u32,
    ) -> (f32, usize) {
        //we calculate all 11, keep only the the max and argmax.
        //we need another one later one, which we calculate by hand.
        //this greatly reduceses memory bandwidth - we needed about 2*11*4 (=88) bytes per position
        //before to keep them all in memory.
        let lls = Coverage::ll(count_a as f32,
            count_c as f32,
            count_g as f32,
            count_t as f32);
        vector_arg_max(&lls)
    }


    fn single_log_likelihood_max_arg_max_plus_other(
        count_a: u32,
        count_c: u32,
        count_g: u32,
        count_t: u32,
        other: usize,
    ) -> ((f32, usize), f32) {
        let lls = Coverage::ll(count_a as f32,
            count_c as f32,
            count_g as f32,
            count_t as f32);
        (vector_arg_max(&lls), lls[other])
    }

    fn score_differences(&self, other: &Self, min_score: f32) -> Vec<ResultRow> {
        let length = self.len();
        if length == 0 || other.len() == 0 {
            return Vec::new();
        }
        let mut result = Vec::new();
        for ii in 0..length {
            // this is a huge speed up for sparse bams.
            // and all rnaseqs are sparse, right
            if (self.0[[ii, BASE_A]] == 0
                && self.0[[ii, BASE_C]] == 0
                && self.0[[ii, BASE_G]] == 0
                && self.0[[ii, BASE_T]] == 0)
                || (other.0[[ii, BASE_A]] == 0
                    && other.0[[ii, BASE_C]] == 0
                    && other.0[[ii, BASE_G]] == 0
                    && other.0[[ii, BASE_T]] == 0) {
                continue;
            }
            let (self_max, self_argmax) = Coverage::single_log_likelihood_max_arg_max(
                self.0[[ii, BASE_A]],
                self.0[[ii, BASE_C]],
                self.0[[ii, BASE_G]],
                self.0[[ii, BASE_T]],
            );
            let ((other_max, other_argmax), other_self_argmax) =
                Coverage::single_log_likelihood_max_arg_max_plus_other(
                    other.0[[ii, BASE_A]],
                    other.0[[ii, BASE_C]],
                    other.0[[ii, BASE_G]],
                    other.0[[ii, BASE_T]],
                    self_argmax,
                );
            if self_argmax == other_argmax {
                // no disagreement
                continue;
            }
            let self_other_argmax = Coverage::single_log_likelihood(
                self.0[[ii, BASE_A]],
                self.0[[ii, BASE_C]],
                self.0[[ii, BASE_G]],
                self.0[[ii, BASE_T]],
                other_argmax,
            );

            let ll_differing = self_max + other_max;
            let ll_same_haplotype_a = self_max + other_self_argmax;
            let ll_same_haplotype_b = self_other_argmax + other_max;
            let ll_same = ll_same_haplotype_a.max(ll_same_haplotype_b);
            let score = ll_differing - ll_same;
            if score >= min_score {
                result.push(ResultRow {
                    relative_pos: ii as u32,
                    count_self_a: self.0[[ii, BASE_A]],
                    count_self_c: self.0[[ii, BASE_C]],
                    count_self_g: self.0[[ii, BASE_G]],
                    count_self_t: self.0[[ii, BASE_T]],
                    count_other_a: other.0[[ii, BASE_A]],
                    count_other_c: other.0[[ii, BASE_C]],
                    count_other_g: other.0[[ii, BASE_G]],
                    count_other_t: other.0[[ii, BASE_T]],
                    haplotype_self: self_argmax as u8,
                    haplotype_other: other_argmax as u8,
                    score,
                });
            }
        }
        result
    }
}

fn is_homo_polymer(seq: &bam::record::Seq, threshold: u8) -> bool {
    let mut counter: u8 = 1;
    let mut last: u8 = 15; //N
    for ii in 0..seq.len() {
        let letter = seq.encoded_base(ii);
        if last == letter {
            counter += 1;
            if counter >= threshold {
                return true;
            }
        } else {
            counter = 1;
        }
        last = letter;
    }
    false
}

fn default_quality_threshold() -> u8 {
    15u8
}

#[derive(Deserialize, Debug)]
struct RunConfig {
    output_dir: String,
    chromosomes: Option<Vec<String>>,
    samples: HashMap<String, Vec<String>>,
    #[serde(default = "default_quality_threshold")]
    quality_threshold: u8,
    filter_homo_polymer_threshold: Option<u8>,
    min_score: Option<f32>
}

fn snp_diff_from_toml(input: &str) -> Result<(), ()> {
    let config: RunConfig = toml::from_str(input).unwrap();
    run_snp_diff(config)
}

fn run_snp_diff(config: RunConfig) -> Result<(), ()> {
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
            .collect();

        let quality_threshold = config.quality_threshold;
        let filter_homo_polymer_threshold = config.filter_homo_polymer_threshold.clone();
        let min_score = config.min_score.unwrap_or(50.0);

        let chunks = chunked_genome::ChunkedGenome::new(first_bam, &config.chromosomes);
        for chunk in chunks.iter(50_000_000) {
            let cov = get_coverages(&config, &chunk, quality_threshold, &filter_homo_polymer_threshold);
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

fn get_coverages(config: &RunConfig, chunk: &chunked_genome::Chunk, quality_threshold: u8, filter_homo_polymer_threshold:&Option<u8>) -> HashMap<String, Coverage> {
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
                .collect()}

fn calculate_differences(coverages: &HashMap<String, Coverage>, pairs: &Vec<Vec<&String>>, output_handles: &mut Vec<std::fs::File>, chunk: &chunked_genome::Chunk, min_score: f32) {
            for (ii, pair) in pairs.iter().enumerate() {
                let delta = coverages[pair[0]].score_differences(&coverages[pair[1]], min_score);
                if !delta.is_empty() {
                    write_results(&mut output_handles[ii], &chunk.chr, chunk.start, delta);
                }
            }
}


fn write_results(
    output_handle: &mut std::fs::File,
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
    for row in rows {
        writeln!(
            output_handle,
            "\t{}\t{}\t{:.13}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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

fn chromosome_name_to_tids(bam: &mut bam::IndexedReader, names: &Vec<String>) -> Vec<u32> {
    let header = bam.header();
    let mut res = Vec::new();
    for name in names {
        res.push(header.tid(name.as_bytes()).unwrap());
    }

    res
}
fn chromosomes_from_bam(bam: &mut bam::IndexedReader) -> Vec<u32> {
    let header = bam.header();
    (0..header.target_count()).collect()
}

fn main() {
    use std::env;
   let args: Vec<String> = env::args().collect();
    let toml = if args.contains(&"--large".to_owned()) {
        //214s as of 10:25 using one core- python needs 818 using
        //183 as of 15:58
"
output_dir = 'tests/ERR329501'
[samples]
    A = ['sample_data/ERR329501.bam']
    B = ['sample_data/GSM1553106.bam']
"
    } else if args.contains(&"--tiny".to_owned())  {
        "
        output_dir = 'tests/test_sample_data'
        [samples]
            A = ['sample_data/sample_a.bam']
            B = ['sample_data/sample_b.bam']
        "

    } else
    { // takes about 50 seconds as of 10:25 // down to 3 seconds at 15:30
"
output_dir = 'tests/ERR329501_chr4'
[samples]
    A = ['sample_data/ERR329501_chr4.bam']
    B = ['sample_data/GSM1553106_chr4.bam']
"
    };
    use std::path::Path;
    let output_dir = Path::new("tests/");
    if output_dir.exists() {
        std::fs::remove_dir_all(output_dir).ok();
    }
    std::fs::create_dir_all(output_dir).unwrap();
    snp_diff_from_toml(toml).unwrap();
    println!("Hello, world!");
}

mod test {

    use super::{Coverage, AA, AC, AG, AT, CC, CG, CT, GG, GT, NN, TT};
    use approx::AbsDiffEq;

    #[test]
    fn test_count_coverage_simple() {
        use super::{Coverage, BASE_A, BASE_C, BASE_G, BASE_T};
        let filename = std::path::Path::new("sample_data/sample_a.bam");
        let start: usize = 10400;
        let cov = Coverage::from_bam(&filename, 0u32, start as u32, (start + 200) as u32, 0, &None);
        assert!(cov.len() == 200);
        assert!(cov.0[[10556 - start, BASE_A]] == 0);
        assert!(cov.0[[10556 - start, BASE_C]] == 0);
        assert!(cov.0[[10556 - start, BASE_G]] == 51);
        assert!(cov.0[[10556 - start, BASE_T]] == 0);
    }
    #[test]
    fn test_homopolymer() {
        use super::is_homo_polymer;
        use rust_htslib::bam;
        let mut read = bam::Record::new();
        read.set(b"A", None, b"AGTC", b"bbbb");
        assert!(!is_homo_polymer(&read.seq(), 5));
        assert!(!is_homo_polymer(&read.seq(), 1)); // 1 makes no sense :)

        let mut read = bam::Record::new();
        read.set(b"A", None, b"AGGTC", b"bbbbb");
        assert!(!is_homo_polymer(&read.seq(), 3));
        assert!(is_homo_polymer(&read.seq(), 2));

        let mut read = bam::Record::new();
        read.set(b"A", None, b"AAGTC", b"bbbbb");
        assert!(!is_homo_polymer(&read.seq(), 3));
        assert!(is_homo_polymer(&read.seq(), 2));

        let mut read = bam::Record::new();
        read.set(b"A", None, b"AAGTTTC", b"AAGTTTC");
        assert!(!is_homo_polymer(&read.seq(), 4));
        assert!(is_homo_polymer(&read.seq(), 3));

        let mut read = bam::Record::new();
        read.set(b"A", None, b"AAGTTTCCCC", b"AAGTTTCCCC");
        assert!(!is_homo_polymer(&read.seq(), 5));
        assert!(is_homo_polymer(&read.seq(), 4));
    }

    #[test]
    fn test_log_likelihood() {
        /*
        use super::{Coverage, HaplotypeLogLikelihoods};
        use ndarray::prelude::*;
        let count_a = vec![100 as u32, 0, 0, 100, 25];
        let count_c = vec![0 as u32, 200, 0, 0, 25];
        let count_g = vec![0 as u32, 0, 100, 0, 25];
        let count_t = vec![0 as u32, 0, 0, 100, 25];

        let mut should = HaplotypeLogLikelihoods::new(5);
        should.0.index_axis_mut(Axis(0), 0).assign(&array![
            -1.00050033e-01,
            -6.94147720e+01,
            -6.94147720e+01,
            -6.94147720e+01,
            -8.00636780e+02,
            -7.60090271e+02,
            -7.60090271e+02,
            -8.00636780e+02,
            -7.60090271e+02,
            -8.00636780e+02,
            -1.38629436e+02
        ]);
        should.0.index_axis_mut(Axis(0), 1).assign(&array![
            -1.60127356e+03,
            -1.38829544e+02,
            -1.52018054e+03,
            -1.52018054e+03,
            -2.00100067e-01,
            -1.38829544e+02,
            -1.38829544e+02,
            -1.60127356e+03,
            -1.52018054e+03,
            -1.60127356e+03,
            -2.77258872e+02
        ]);
        should.0.index_axis_mut(Axis(0), 2).assign(&array![
            -8.00636780e+02,
            -7.60090271e+02,
            -6.94147720e+01,
            -7.60090271e+02,
            -8.00636780e+02,
            -6.94147720e+01,
            -7.60090271e+02,
            -1.00050033e-01,
            -6.94147720e+01,
            -8.00636780e+02,
            -1.38629436e+02
        ]);
        should.0.index_axis_mut(Axis(0), 3).assign(&array![
            -8.00736830e+02,
            -8.29505066e+02,
            -8.29505066e+02,
            -1.38829544e+02,
            -1.60127356e+03,
            -1.52018054e+03,
            -8.29505066e+02,
            -1.60127356e+03,
            -8.29505066e+02,
            -8.00736830e+02,
            -2.77258872e+02
        ]);
        should.0.index_axis_mut(Axis(0), 4).assign(&array![
            -6.00502597e+02,
            -4.14752502e+02,
            -4.14752533e+02,
            -4.14752533e+02,
            -6.00502597e+02,
            -4.14752533e+02,
            -4.14752533e+02,
            -6.00502597e+02,
            -4.14752533e+02,
            -6.00502597e+02,
            -1.38629436e+02
        ]);
        let cov = Coverage::from_counts(count_a, count_c, count_g, count_t);
        let actual = cov.log_likelihood();
        assert!(actual.0.abs_diff_eq(&should.0, 1e-4))
           */
    }

    #[test]
    fn test_score_differences() {
        let cov_a = Coverage::from_counts(
            vec![100, 0, 0, 25],
            vec![0, 100, 0, 25],
            vec![0, 0, 100, 25],
            vec![0, 0, 100, 25],
        );
        let cov_b = Coverage::from_counts(
            vec![0, 0, 100, 50],
            vec![100, 0, 0, 50],
            vec![0, 0, 0, 0],
            vec![0, 0, 100, 0],
        );
        let res = cov_a.score_differences(&cov_b, 0.0);
        assert!(res.len() == 3);
        assert_eq!(res[0].relative_pos, 0);
        assert_eq!(res[1].relative_pos, 2);
        assert_eq!(res[2].relative_pos, 3);
        assert!(res[0].score.abs_diff_eq(&800.53764f32, 1e-3));
        assert!(res[1].score.abs_diff_eq(&690.67554f32, 1e-3));
        assert!(res[2].score.abs_diff_eq(&69.21466f32, 1e-3));
        assert_eq!(res[0].haplotype_self, AA as u8);
        assert_eq!(res[0].haplotype_other, CC as u8);
        assert_eq!(res[1].haplotype_self, GT as u8);
        assert_eq!(res[1].haplotype_other, AT as u8);
        assert_eq!(res[2].haplotype_self, NN as u8);
        assert_eq!(res[2].haplotype_other, AC as u8);
    }

    #[test]
    fn test_score_differences_problematic() {
        use super::Coverage;
        use approx::AbsDiffEq;
        let cov_a = Coverage::from_counts(vec![0, 51], vec![0, 0], vec![51, 0], vec![0, 0]);
        let cov_b = Coverage::from_counts(vec![0, 0], vec![0, 0], vec![0, 0], vec![51, 51]);
        let res = cov_a.score_differences(&cov_b, 0.);
        assert_eq!(res[0].score, res[1].score);
        assert!(res[0].score.abs_diff_eq(&408.2737f32, 1e-4))
    }

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
