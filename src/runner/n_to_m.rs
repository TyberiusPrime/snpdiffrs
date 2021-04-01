use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use progressing::{mapping::Bar as MappingBar, Baring};

#[allow(unused_imports)]
use slog::{debug, info};

use super::{all_files_exists, ensure_output_dir, get_logger, haplotype_const_to_str, RunConfig};
use super::runner::{Runner, RunStrategy};



struct NtoMStrategy {
    queries: Vec<String>, //something set-ish?
    references: Vec<String>,
}

impl NtoMStrategy {
    fn new(config: &RunConfig) -> Self {
        NtoMStrategy{
            queries: config.queries.as_ref().expect("No queries set").keys().map(|x|x.to_string()).collect(),
            references: config.references.as_ref().expect("No queries set").keys().map(|x|x.to_string()).collect(),
        }
    }

}

impl RunStrategy for NtoMStrategy {
    fn expected_number_of_comparisons(&self) -> usize {
        self.queries.len() * self.references.len()
    }

    fn expected_number_of_uses_per_block(&self, sample_name: &str) -> usize {
        if self.queries.iter().any(|x| x == sample_name) {
            self.references.len()
        }
        else {
            self.queries.len()
        }
    }

    fn include_comparison(&self, sample_a_name: &str, sample_b_name: &str) -> bool {
        (self.queries.iter().any(|x| x == sample_a_name)  && self.references.iter().any(|x| x == sample_b_name)) ||
        (self.queries.iter().any(|x| x == sample_b_name)  && self.references.iter().any(|x| x == sample_a_name))
    }
}

pub fn run_n_to_m(config: RunConfig) -> Result<(), ()> {
    all_files_exists(&config.queries, "queries");
    all_files_exists(&config.references, "references");
    ensure_output_dir(&Path::new(&config.output_dir));

    let samples_and_input_filenames = config.queries_references_and_input_filenames();
    let log = get_logger();
    debug!(log, "Samples&input_files: {:?}", samples_and_input_filenames);
    let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) =
            samples_and_input_filenames.into_iter().unzip();
    //todo: ensure they are distinct!
    let strat = Box::new(NtoMStrategy::new(&config));
    let runner = Runner::new(config,
        samples,
        input_filenames,
        strat);
    runner.run();
    Ok(())
}


