use std::path::Path;
use std::convert::TryInto;
use rust_htslib::bam;
use rust_htslib::htslib;
//use rust_htslib::prelude::*;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;


fn main() {
    println!("Hello, world!");
}

struct Coverage{
    a: Vec<u32>,
    c: Vec<u32>,
    g: Vec<u32>,
    t: Vec<u32>,
}

impl Coverage {
    fn new(length: usize) -> Coverage{
        Coverage{
        a: vec![0; length],
        c: vec![0; length],
        g: vec![0; length],
        t: vec![0; length],

    }
    }
}


fn count_coverage<P: AsRef<Path>>(filename: P, chr: &str, start: i64, stop: i64, quality_threshold: u8) 
-> Coverage {
    let mut bam = bam::IndexedReader::from_path(filename).expect("Could not read input bam");
    let tid: u32 = bam.header().tid(chr.as_bytes()).expect("Could not find chromosome");
    dbg!((tid, start, stop));
    bam.fetch(tid, start as u64, stop as u64).unwrap();
    let length: u64 = (stop - start).try_into().expect("stop < start");
    dbg!(length);
    let mut result: Coverage = Coverage::new(length as usize);
    let mut read: bam::Record = bam::Record::new();
    while let Ok(true) = bam.read(&mut read) {
        if (read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY | htslib::BAM_FQCFAIL | htslib::BAM_FDUP) as u16) > 0{
         continue
        }
        let seq = read.seq();
        for [read_pos, genome_pos] in read.aligned_pairs().iter() {
            if (*genome_pos < start) || *genome_pos >= stop{
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
            if !read.is_reverse() {
                match base {
                    1 => result.a[(genome_pos - start) as usize] += 1,
                    2 => result.c[(genome_pos - start) as usize] += 1,
                    4 => result.g[(genome_pos - start) as usize] += 1,
                    8 => result.t[(genome_pos - start) as usize] += 1,
                    _ => {}
                };
            } else {
                match base {
                    2 => result.a[(genome_pos - start) as usize] += 1,
                    1 => result.c[(genome_pos - start) as usize] += 1,
                    8 => result.g[(genome_pos - start) as usize] += 1,
                    4 => result.t[(genome_pos - start) as usize] += 1,
                    _ => {}
                };
            }

        }
    }
    result
}

mod test {
    use std::path::Path;
    #[test]
    fn test_count_coverage_simple() {
        use super::count_coverage;
        let filename = Path::new("sample_data/sample_a.bam");
        let start: usize = 10400;
        let cov = count_coverage(&filename, "1", start as i64,(start + 200) as i64, 0);
        dbg!(&cov.a[10556 - start]);
        dbg!(&cov.c[10556 - start]);
        dbg!(&cov.g[10556 - start]);
        dbg!(&cov.t[10556 - start]);
        assert!(cov.a.len() == 200);
        assert!(cov.a[10556 - start] == 0);
        assert!(cov.c[10556 - start] == 0);
        assert!(cov.g[10556 - start] == 0);
        assert!(cov.t[10556 - start] == 51);
    }
}
