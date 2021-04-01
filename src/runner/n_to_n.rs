
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


struct NtoNStrategy {
    no_of_samples: usize
}

impl NtoNStrategy {
    fn new(config: &RunConfig) -> Self {
        NtoNStrategy{no_of_samples: config.samples.as_ref().expect("No samples provided").len()}
    }

}

impl RunStrategy for NtoNStrategy {
    fn expected_number_of_comparisons(&self) -> usize {
        self.no_of_samples * (self.no_of_samples -1) / 2
    }

    fn expected_number_of_uses_per_block(&self, sample_name: &str) -> usize {
        self.no_of_samples - 1 // against n-1 other blocks, but not against it self
    }

    fn include_comparison(&self, sample_a_name: &str, sample_b_name: &str) -> bool {
        //we want all the pairs, but only one of them
        use std::cmp::Ordering;
        match sample_a_name.cmp(sample_b_name) {
            Ordering::Less => true,
            _ => false
        }
    }
}

pub fn run_n_to_n(config: RunConfig) -> Result<(), ()> {
    all_files_exists(&config.samples, "samples");
    ensure_output_dir(&Path::new(&config.output_dir));

    let samples_and_input_filenames = config.samples_and_input_filenames();
    let (samples, input_filenames): (Vec<String>, Vec<Vec<PathBuf>>) =
            samples_and_input_filenames.into_iter().unzip();
    let strat = Box::new(NtoNStrategy::new(&config));
    let runner = Runner::new(config,
        samples,
        input_filenames,
        strat);
    runner.run();
    Ok(())
}

mod test {
    #[cfg(test)]
    use super::super::snp_diff_from_toml;

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
        snp_diff_from_toml(toml).unwrap();
        let should =
            std::fs::read_to_string("sample_data/marsnpdiff_sample_a_vs_sample_b.tsv").unwrap();
        let should = should.replace(".0", "");
        let actual = std::fs::read_to_string("tests/test_sample_data/A_vs_B.tsv").unwrap();
        assert_eq!(should, actual)
    }
    #[test]
    #[should_panic]
    fn test_same_data() {
        let toml = "
output_dir = 'tests/test_same_data'
[samples]
    A = ['sample_data/sample_a.bam']
    A = ['sample_data/sample_a.bam']
";
        use std::path::Path;
        let output_dir = Path::new("tests/test_same_data");
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        snp_diff_from_toml(toml).unwrap();
    }

    #[test]
    fn test_sample_n_to_n() {
        let test_path = "tests/test_sample_data_n_to_n";
        let toml = format!(
            "
output_dir = '{}'
[samples]
    A = ['sample_data/sample_a.bam'] # mixed up order is intentional
    C = ['sample_data/sample_a.bam']
    D = ['sample_data/sample_b.bam']
    B = ['sample_data/sample_b.bam']
",
            test_path
        );
        use std::path::Path;
        let output_dir = Path::new(test_path);
        if output_dir.exists() {
            std::fs::remove_dir_all(output_dir).ok();
        }
        std::fs::create_dir_all(output_dir).unwrap();
        snp_diff_from_toml(&toml).unwrap();
        let should =
            std::fs::read_to_string("sample_data/marsnpdiff_sample_a_vs_sample_b.tsv").unwrap();
        let should = should.replace(".0", "");

        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_B.tsv").unwrap();
        assert_eq!(should, actual);
        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_D.tsv").unwrap();
        assert_eq!(should, actual);

        let actual = std::fs::read_to_string(test_path.to_owned() + "/A_vs_C.tsv").unwrap();
        assert_eq!(actual.matches("\n").count(), 1);

        let actual = std::fs::read_to_string(test_path.to_owned() + "/B_vs_D.tsv").unwrap();
        assert_eq!(actual.matches("\n").count(), 1);
        std::fs::remove_dir_all(output_dir).ok();
    }
}
