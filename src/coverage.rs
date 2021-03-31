use crate::consts::*;
use lzzzz::{lz4, lz4_hc, lz4f};
use ndarray::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::htslib;
use serde::{Deserialize, Serialize};
use std::convert::TryInto;
use std::path::Path;

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

fn vector_arg_max(input: &[f32; 11]) -> (f32, usize) {
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
/// Central 'count and store Coverage type
/// we only count to u16 (65k) to same on memory.
/// At those sequencing depths we overrun our f32
/// probably anyhow.
pub struct Coverage(pub Array2<u16>);

impl std::fmt::Debug for Coverage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Coverage")
    }
}

#[derive(Debug)]
pub struct ResultRow {
    pub relative_pos: u32,
    pub count_self_a: u16,
    pub count_self_c: u16,
    pub count_self_g: u16,
    pub count_self_t: u16,
    pub count_other_a: u16,
    pub count_other_c: u16,
    pub count_other_g: u16,
    pub count_other_t: u16,
    pub score: f32,
    pub haplotype_self: u8,
    pub haplotype_other: u8,
}

fn u32_to_u16_saturating(input: u32) -> u16 {
    if input > u16::MAX as u32 {
        u16::MAX
    } else {
        input as u16
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct EncCoverageEntry {
    pub offset: u32,
    pub count_a: u16,
    pub count_c: u16,
    pub count_g: u16,
    pub count_t: u16,
}

#[derive(Serialize, Deserialize)]
pub struct EncCoverage{
    pub entries: Vec<EncCoverageEntry>,
    length: usize,

}

impl EncCoverage {
    pub fn new(input: &Coverage) -> EncCoverage {
        let mut res = Vec::new();
        for (ii, row) in input.0.axis_iter(Axis(0)).enumerate()
        {
            if row[0] != 0 ||
                row[1] != 0 ||
                row[2] != 0 ||
                row[3] != 0 {
                    res.push(
                        EncCoverageEntry{
                            offset: ii as u32,
                            count_a: row[0],
                            count_c: row[1],
                            count_g: row[2],
                            count_t: row[3],
                });
            }
        }
        EncCoverage{entries: res, length: input.0.len()}
    }

    pub fn from_preprocessed(preprocessed_file: &Path) -> Option<Self> {
        use byteorder::ByteOrder;
        use byteorder::{LittleEndian, ReadBytesExt};

        let fh = std::fs::File::open(preprocessed_file).expect("Could not open file");
        let mut fh = std::io::BufReader::new(fh);
        let mut raw = Vec::new();
        zstd::stream::copy_decode(&mut fh, &mut raw).expect("decompressing failed");
        let mut raw = &raw[..];

        let mut result: Vec<EncCoverageEntry> = Vec::new();
        let length = raw.read_u64::<LittleEndian>().unwrap();
        loop {
            let offset = match raw.read_u32::<LittleEndian>() {
                Ok(o) => o,
                Err(_) => {break;}
            };
            result.push(
                EncCoverageEntry {
                    offset: offset,
                    count_a: raw.read_u16::<LittleEndian>().unwrap(),
                    count_c: raw.read_u16::<LittleEndian>().unwrap(),
                    count_g: raw.read_u16::<LittleEndian>().unwrap(),
                    count_t: raw.read_u16::<LittleEndian>().unwrap(),
                });
        };
        Some(EncCoverage{
            entries: result,
            length: length as usize,
        })
    }

    pub fn sum(&self) -> usize {
        let mut total: usize = 0;
        for row in self.entries.iter() {
            total += row.count_a as usize;
            total += row.count_c as usize;
            total += row.count_g as usize;
            total += row.count_t as usize;
        }
        total
    }

    pub fn len(&self)  -> usize{
        self.length
    }

    fn score_single_difference(mine: &EncCoverageEntry, other: &EncCoverageEntry, result: &mut Vec<ResultRow>, min_score: f32, offset: u32) {
        // this is a huge speed up for sparse bams.
        // and all rnaseqs are sparse, right
        let sa = mine.count_a;
        let sc = mine.count_c;
        let sg = mine.count_g;
        let st = mine.count_t;
        if sa == 0 && sc == 0 && sg == 0 && st == 0 {
            panic!("should not happen");
        }
        let oa = other.count_a;
        let oc = other.count_c;
        let og = other.count_g;
        let ot = other.count_t;
        if oa == 0 && oc == 0 && og == 0 && ot == 0 {
            panic!("should not happen2");
        }
        let (self_max, self_argmax) =
            Coverage::single_log_likelihood_max_arg_max(sa, sc, sg, st);
        let ((other_max, other_argmax), other_self_argmax) =
            Coverage::single_log_likelihood_max_arg_max_plus_other(oa, oc, og, ot, self_argmax);
        if self_argmax == other_argmax {
            // no disagreement
            return;
        }
        let self_other_argmax = Coverage::single_log_likelihood(sa, sc, sg, st, other_argmax);

        let ll_differing = self_max + other_max;
        let ll_same_haplotype_a = self_max + other_self_argmax;
        let ll_same_haplotype_b = self_other_argmax + other_max;
        let ll_same = ll_same_haplotype_a.max(ll_same_haplotype_b);
        let score = ll_differing - ll_same;
        if score >= min_score {
            result.push(ResultRow {
                relative_pos: mine.offset as u32 + offset,
                count_self_a: sa,
                count_self_c: sc,
                count_self_g: sg,
                count_self_t: st,
                count_other_a: oa,
                count_other_c: oc,
                count_other_g: og,
                count_other_t: ot,
                haplotype_self: self_argmax as u8,
                haplotype_other: other_argmax as u8,
                score,
            });
        }
    }

    pub fn score_differences(&self, other: &Self, min_score: f32, offset: u32) -> Vec<ResultRow> {
        use std::cmp::Ordering;
        if self.len() == 0 || other.len() == 0 {
            return Vec::new();
        }
        let mut iter_mine = self.entries.iter();
        let mut iter_other = other.entries.iter();
        let mut element_mine = iter_mine.next();
        let mut element_other = iter_other.next();
        let mut res = Vec::new();
        loop {
            match (element_mine, element_other) {
                (None, None) => return res,
                (None, _) => return res,
                (_, None,)  => return res,
                (Some(m), Some(o)) => {
                    match m.offset.cmp(&o.offset) {
                        Ordering::Equal => {
                            EncCoverage::score_single_difference(m, o, &mut res, min_score, offset);
                            element_mine = iter_mine.next();
                            element_other = iter_other.next();
                        },
                        Ordering::Less => element_mine = iter_mine.next(),
                        Ordering::Greater => element_other = iter_other.next(),
                        }
                    }

                }
            }
    }
}

impl std::fmt::Debug for EncCoverage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("EncCoverage")
    }
}



impl Coverage {
    pub fn new(length: usize) -> Coverage {
        Coverage(Array2::zeros((length, 4)))
    }

    /*pub fn clear(&mut self) {
        self.0.fill(0);
    }
    */

    #[allow(dead_code)]
    pub fn from_counts(
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
            res[[ii, 0]] = u32_to_u16_saturating(count_a[ii]);
            res[[ii, 1]] = u32_to_u16_saturating(count_c[ii]);
            res[[ii, 2]] = u32_to_u16_saturating(count_g[ii]);
            res[[ii, 3]] = u32_to_u16_saturating(count_t[ii]);
        }
        Coverage(res)
    }

    #[allow(dead_code)]
    pub fn from_bam(
        filename: &Path,
        tid: u32,
        start: u32,
        stop: u32,
        quality_threshold: u8,
        filter_homo_polymer_threshold: &Option<u8>,
    ) -> Coverage {
        Coverage::from_bams(
            &[filename],
            tid,
            start,
            stop,
            quality_threshold,
            filter_homo_polymer_threshold,
        )
    }

    pub fn from_bams(
        bams: &[&Path],
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
            // whether we have any reads at all.
            result
        } else {
            Coverage::new(0)
        }
    }

    pub fn sum(&self) -> usize {
        let mut total: usize = 0;
        for i in 0..self.0.shape()[0] {
            total += self.0[(i, 0)] as usize;
            total += self.0[(i, 1)] as usize;
            total += self.0[(i, 2)] as usize;
            total += self.0[(i, 3)] as usize;
        }
        total
    }

    pub fn head(&self) {
        use ndarray::Slice;

        let b = self.0.slice_axis(Axis(0), Slice::from(..100));
        println!("{:?}", b);
    }

    pub fn first_delta(&self, other: &Self) {
        use ndarray::Slice;
        for row in 0..self.0.shape()[0] {
            let a = self.0.slice_axis(Axis(0), Slice::from(row..row + 1));
            let b = other.0.slice_axis(Axis(0), Slice::from(row..row + 1));
            if a != b {
                println!("Difference row {},\n {:?}\n{:?}", row, a, b);
                break;
            }
        }
    }

    //read start,stop, genome start, stop aligned blocks
    //replaces aligned_pairs with just the ends
    //and such safes on allocations.
    fn aligned_blocks(record: &bam::Record) -> Vec<(i64, i64, i64, i64)> {
        use rust_htslib::bam::record::Cigar;
        let mut result = Vec::new();

        let mut pos: i64 = record.pos();
        let mut qpos: i64 = 0;
        for entry in record.cigar().iter() {
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let qstart = qpos;
                    let qend = qstart + *len as i64;
                    let rstart = pos;
                    let rend = pos + *len as i64;
                    result.push((qstart, qend, rstart, rend));
                    qpos += *len as i64;
                    pos += *len as i64;
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    qpos += *len as i64;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    pos += *len as i64;
                }
                Cigar::HardClip(_) => {} // no advance
                Cigar::Pad(_) => panic!("Padding (Cigar::Pad) is not supported."), //padding is only used for multiple sequence alignment
            }
        }
        result
    }

    fn update_from_bam<P: AsRef<Path> + std::fmt::Debug + Clone>(
        &mut self,
        filename: P,
        tid: u32,
        start: u32,
        stop: u32,
        quality_threshold: u8,
        filter_homo_polymer_threshold: &Option<u8>,
    ) -> bool {
        let mut bam = bam::IndexedReader::from_path(filename).expect("Could not read input bam");
        let start = start as i64;
        let stop = stop as i64;
        bam.fetch((tid, start as u64, stop as u64)).unwrap();
        let mut any = false;
        for read in bam.rc_records().filter_map(|x| x.ok()).filter(|read| {
            (read.flags()
                & (htslib::BAM_FUNMAP // 0x4
                    | htslib::BAM_FSECONDARY // 256 = 0x100
                    | htslib::BAM_FQCFAIL // 521 = 0x200
                    | htslib::BAM_FDUP) as u16) // 0x400
                == 0
        }) {
            let seq = read.seq();
            if filter_homo_polymer_threshold.is_some()
                && is_homo_polymer(&seq, filter_homo_polymer_threshold.unwrap())
            {
                continue;
            }
            //for [read_pos, genome_pos] in read.aligned_pairs().iter() {
            for (qstart, qend, genome_start, genome_end) in Coverage::aligned_blocks(&read) {
                if (genome_end < start) || genome_start >= stop {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                for (read_pos, genome_pos) in (qstart..qend).zip(genome_start..genome_end) {
                    if (genome_pos - start) as usize > self.len() - 1 {
                        continue;
                    }
                    //if read_pos as usize > read.qual().len() - 1 {
                    //continue;
                    //}
                    if read.qual()[read_pos as usize] <= quality_threshold {
                        continue;
                    }
                    any = true;
                    let base: u8 = seq.encoded_base((read_pos) as usize);
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
                    self.0[[(genome_pos - start) as usize, out_base]] =
                        self.0[[(genome_pos - start) as usize, out_base]].saturating_add(1);
                }
            }
        }
        any
    }

    pub fn len(self: &Self) -> usize {
        self.0.dim().0
    }

    pub fn single_log_likelihood(
        count_a: u16,
        count_c: u16,
        count_g: u16,
        count_t: u16,
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

    fn ll(a: f32, c: f32, g: f32, t: f32) -> [f32; 11] {
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
            a_ll_495 + c_ll_495 + g_ll_005 + t_ll_005,   //AC
            a_ll_495 + c_ll_005 + g_ll_495 + t_ll_005,   //AG
            a_ll_495 + c_ll_005 + g_ll_005 + t_ll_495,   //AT
            a_ll_005 + c_ll_495 + g_ll_495 + t_ll_005,   //CG
            a_ll_005 + c_ll_495 + g_ll_005 + t_ll_495,   //CT
            a_ll_005 + c_ll_005 + g_ll_495 + t_ll_495,   //GT
            a * *LL_25 + c * *LL_25 + g * *LL_25 + t * *LL_25,
        ]
    }

    fn single_log_likelihood_max_arg_max(
        count_a: u16,
        count_c: u16,
        count_g: u16,
        count_t: u16,
    ) -> (f32, usize) {
        //we calculate all 11, keep only the the max and argmax.
        //we need another one later one, which we calculate by hand.
        //this greatly reduceses memory bandwidth - we needed about 2*11*4 (=88) bytes per position
        //before to keep them all in memory.
        let lls = Coverage::ll(
            count_a as f32,
            count_c as f32,
            count_g as f32,
            count_t as f32,
        );
        vector_arg_max(&lls)
    }

    fn single_log_likelihood_max_arg_max_plus_other(
        count_a: u16,
        count_c: u16,
        count_g: u16,
        count_t: u16,
        other: usize,
    ) -> ((f32, usize), f32) {
        let lls = Coverage::ll(
            count_a as f32,
            count_c as f32,
            count_g as f32,
            count_t as f32,
        );
        (vector_arg_max(&lls), lls[other])
    }

    pub fn score_differences(&self, other: &Self, min_score: f32, offset: u32) -> Vec<ResultRow> {
        let length = self.len();
        if length == 0 || other.len() == 0 {
            return Vec::new();
        }
        let mut result = Vec::new();
        for (ii, (self_row, other_row)) in self
            .0
            .axis_iter(Axis(0))
            .zip(other.0.axis_iter(Axis(0)))
            .enumerate()
        {
            // this is a huge speed up for sparse bams.
            // and all rnaseqs are sparse, right
            let sa = self_row[BASE_A];
            let sc = self_row[BASE_C];
            let sg = self_row[BASE_G];
            let st = self_row[BASE_T];
            if sa == 0 && sc == 0 && sg == 0 && st == 0 {
                continue;
            }
            let oa = other_row[BASE_A];
            let oc = other_row[BASE_C];
            let og = other_row[BASE_G];
            let ot = other_row[BASE_T];
            if oa == 0 && oc == 0 && og == 0 && ot == 0 {
                continue;
            }
            let (self_max, self_argmax) =
                Coverage::single_log_likelihood_max_arg_max(sa, sc, sg, st);
            let ((other_max, other_argmax), other_self_argmax) =
                Coverage::single_log_likelihood_max_arg_max_plus_other(oa, oc, og, ot, self_argmax);
            if self_argmax == other_argmax {
                // no disagreement
                continue;
            }
            let self_other_argmax = Coverage::single_log_likelihood(sa, sc, sg, st, other_argmax);

            let ll_differing = self_max + other_max;
            let ll_same_haplotype_a = self_max + other_self_argmax;
            let ll_same_haplotype_b = self_other_argmax + other_max;
            let ll_same = ll_same_haplotype_a.max(ll_same_haplotype_b);
            let score = ll_differing - ll_same;
            if score >= min_score {
                result.push(ResultRow {
                    relative_pos: ii as u32 + offset,
                    count_self_a: sa,
                    count_self_c: sc,
                    count_self_g: sg,
                    count_self_t: st,
                    count_other_a: oa,
                    count_other_c: oc,
                    count_other_g: og,
                    count_other_t: ot,
                    haplotype_self: self_argmax as u8,
                    haplotype_other: other_argmax as u8,
                    score,
                });
            }
        }
        result
    }

    pub fn into_raw_vec(self) -> Vec<u16> {
        self.0.into_raw_vec() //todo: this is broken by drop
    }

    pub fn equal(&self, other: &Self) -> bool {
        (self.0 == other.0) && (self.len() == other.len())
    }

    pub fn shape(&self) -> &[usize] {
        self.0.shape()
    }
}

/*
impl Drop for Coverage {
    fn drop(&mut self) {
        use slog::debug;
        use crate::runner::get_logger;
        let log = get_logger();
        debug!(log, "Dropped coverage");
    }
}
*/

pub(crate) struct BaseCounts {
    pub a: u64,
    pub g: u64,
    pub c: u64,
    pub t: u64,
}

impl BaseCounts {
    pub(crate) fn new_histogram() -> Vec<BaseCounts> {
        let mut res = Vec::new();
        for __ii in 0..u8::MAX {
            res.push(BaseCounts {
                a: 0,
                c: 0,
                g: 0,
                t: 0,
            });
        }
        res
    }

    pub(crate) fn update(this: &mut Vec<BaseCounts>, other: &[BaseCounts]) {
        for ii in 0..this.len() {
            this[ii].a += other[ii].a;
            this[ii].c += other[ii].c;
            this[ii].g += other[ii].g;
            this[ii].t += other[ii].t;
        }
    }
}

pub(crate) fn count_homopolymers_in_bam(
    filename: &Path,
    tid: u32,
    start: u32,
    stop: u32,
) -> Vec<BaseCounts> {
    let mut res = BaseCounts::new_histogram();
    let mut bam = bam::IndexedReader::from_path(filename).expect("Could not read input bam");
    let start = start as i64;
    let stop = stop as i64;
    bam.fetch((tid, start as u64, stop as u64)).unwrap();
    for read in bam.rc_records().filter_map(|x| x.ok()).filter(|read| {
        (read.flags()
                & (htslib::BAM_FUNMAP // 0x4
                    | htslib::BAM_FSECONDARY // 256 = 0x100
                    | htslib::BAM_FQCFAIL // 521 = 0x200
                    | htslib::BAM_FDUP) as u16) // 0x400
                == 0
    }) {
        if read.pos() < start || read.pos() >= stop {
            continue;
        }
        let (count, base) = quantify_homo_polymer(&read.seq());
        let count = count as usize;
        match base {
            BASE_A => res[count].a += 1,
            BASE_C => res[count].c += 1,
            BASE_G => res[count].g += 1,
            BASE_T => res[count].t += 1,
            _ => panic!("None AGTC homopolymer?"), //ignore non regular base - must be one weird read?
        }
    }
    res
}

fn quantify_homo_polymer(seq: &bam::record::Seq) -> (u8, usize) {
    let mut max_value: u8 = 1;
    let mut max_base = BASE_A;
    let mut current_value: u8 = 0;
    let mut last: u8 = 15; //N
    for ii in 0..seq.len() {
        let letter = seq.encoded_base(ii);
        if last == letter {
            current_value = current_value.saturating_add(1);
            if current_value >= max_value {
                max_value = current_value;
                max_base = match letter {
                    1 => BASE_A,
                    2 => BASE_C,
                    4 => BASE_G,
                    8 => BASE_T,
                    _ => continue,
                };
            }
        } else {
            current_value = 1;
        }
        last = letter;
    }
    (max_value, max_base)
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

mod test {

    #[allow(unused_imports)]
    use super::{Coverage, AA, AC, AG, AT, CC, CG, CT, GG, GT, NN, TT};
    #[cfg(test)]
    use approx::AbsDiffEq;

    #[test]
    fn test_count_coverage_simple() {
        use super::{Coverage, BASE_A, BASE_C, BASE_G, BASE_T};
        let filename = std::path::Path::new("sample_data/sample_a.bam");
        let start: usize = 10400;
        let cov = Coverage::from_bam(
            &filename,
            0u32,
            start as u32,
            (start + 200) as u32,
            0,
            &None,
        );
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
        let res = cov_a.score_differences(&cov_b, 0.0, 0);
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
        let cov_a = Coverage::from_counts(vec![0, 51], vec![0, 0], vec![51, 0], vec![0, 0]);
        let cov_b = Coverage::from_counts(vec![0, 0], vec![0, 0], vec![0, 0], vec![51, 51]);
        let res = cov_a.score_differences(&cov_b, 0., 0);
        assert_eq!(res[0].score, res[1].score);
        assert!(res[0].score.abs_diff_eq(&408.2737f32, 1e-4))
    }
}
