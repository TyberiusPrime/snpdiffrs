use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::str;

pub struct ChunkedGenome {
    bam: bam::IndexedReader,
    chromosomes: Vec<String>,
}

impl ChunkedGenome {
    ///create a new chunked genome for iteration
    ///if you pass in a tree, it is guranteed that the splits happen
    ///between entries of the tree, not inside.
    pub fn new(bam: bam::IndexedReader, chromosomes: &Option<Vec<String>>) -> ChunkedGenome {
        ChunkedGenome {
            chromosomes: match chromosomes {
                Some(x) => {
                    let chromosomes: Vec<String> = x.iter().map(|y| y.to_string()).collect();
                    let targets = bam.header().target_names();
                    for c in chromosomes.iter() {
                        if !targets.contains(&c.as_bytes()) {
                            panic!("invalid chromosome specified: {}", c);
                        }
                    }
                    chromosomes
                }
                None => bam
                    .header()
                    .target_names()
                    .iter()
                    .map(|x| str::from_utf8(x).unwrap().to_string())
                    .collect(),
            },
            bam,
        }
    }
    pub fn iter(&self, chunk_size: usize) -> ChunkedGenomeIterator {
        ChunkedGenomeIterator {
            cg: &self,
            it: self.chromosomes.iter(),
            last_start: 0,
            last_tid: 0,
            last_chr_length: 0,
            last_chr: "".to_string(),
            chunk_size: chunk_size as u32,
        }
    }
}

pub struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome,
    it: std::slice::Iter<'a, String>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u64,
    chunk_size: u32,
}
#[derive(Debug, Clone)]
pub struct Chunk {
    pub chr: String,
    pub tid: u32,
    pub start: u32,
    pub stop: u32,
}

impl<'a> Iterator for ChunkedGenomeIterator<'a> {
    type Item = Chunk;
    fn next(&mut self) -> Option<Chunk> {
        let chunk_size = self.chunk_size;
        if self.last_start as u64 >= self.last_chr_length {
            let next_chr = match self.it.next() {
                Some(x) => x,
                None => return None,
            };
            let tid = self.cg.bam.header().tid(next_chr.as_bytes()).unwrap();
            let chr_length = self.cg.bam.header().target_len(tid).unwrap();
            self.last_tid = tid;
            self.last_chr_length = chr_length;
            self.last_chr = next_chr.to_string();
            self.last_start = 0;
        }

        let stop = self.last_start + chunk_size;
        let c = Chunk {
            chr: self.last_chr.clone(),
            tid: self.last_tid,
            start: self.last_start,
            stop,
        };
        self.last_start = stop;
        Some(c)
    }
}
