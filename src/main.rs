mod chunked_genome;
mod consts;
mod coverage;
mod runner;

#[macro_use]
extern crate lazy_static;

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
    } else if args.contains(&"--tiny".to_owned()) {
        "
        output_dir = 'tests/test_sample_data'
        [samples]
            A = ['sample_data/sample_a.bam']
            B = ['sample_data/sample_b.bam']
        "
    }  else if args.contains(&"--mt".to_owned()) {
        "
        output_dir = 'tests/test_sample_data'
        chromosomes = ['MT']
        min_score = 0.0
        [samples]
            A = ['sample_data/ERR329501.bam']
            B = ['sample_data/GSM1553106.bam']
        "
    } else {
        // takes about 50 seconds as of 10:25 // down to 3 seconds at 15:30
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
    runner::snp_diff_from_toml(toml).unwrap();
    println!("Hello, world!");
}
