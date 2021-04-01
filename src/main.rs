//todo: remove these, clean up
#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
mod chunked_genome;
mod consts;
mod coverage;
mod runner;

#[macro_use]
extern crate lazy_static;
use libc::{prctl, PR_SET_PDEATHSIG};
use nix::sys::signal::SIGTERM;
use std::panic;

fn main() {
    //die if the parent process is killed.
    unsafe {
        prctl(PR_SET_PDEATHSIG, SIGTERM);
    }

    //die if a panic occurs in any thread
    let orig_hook = panic::take_hook();
    panic::set_hook(Box::new(move |panic_info| {
        // invoke the default handler and exit the process
        orig_hook(panic_info);
        std::process::exit(1);
    }));

    use std::env;
    let args: Vec<String> = env::args().collect();
    let toml = if args.contains(&"--large".to_owned()) {
        //todo: replace with examples
        //214s as of 10:25 using one core- python needs 818 using
        //183 as of 15:58
        //down to 25s with multi core and latest optimizations as of 08-09-2020
        "
output_dir = 'tests/ERR329501'
[samples]
    A = ['sample_data/ERR329501.bam']
    B = ['sample_data/GSM1553106.bam']
"
        .to_string()
    } else if args.contains(&"--tiny".to_owned()) {
        "
        output_dir = 'tests/test_sample_data'
        [samples]
            A = ['sample_data/sample_a.bam']
            B = ['sample_data/sample_b.bam']
        "
        .to_string()
    } else if args.contains(&"--mt".to_owned()) {
        "
        output_dir = 'tests/test_sample_data'
        chromosomes = ['MT']
        min_score = 0.0
        [samples]
            A = ['sample_data/ERR329501.bam']
            B = ['sample_data/GSM1553106.bam']
        "
        .to_string()
    } else if args.contains(&"--normal".to_owned()) {
        // takes about 50 seconds as of 10:25 // down to 3 seconds at 15:30
        "
            output_dir = 'tests/ERR329501_chr4'
            [samples]
                A = ['sample_data/ERR329501_chr4.bam']
                B = ['sample_data/GSM1553106_chr4.bam']
            "
        .to_string()
    } else {
        if args.len() < 2 {
            panic!("Please pass in a configuration file"); //todo: pretty help
        }
        let filename = &args[1];
        std::fs::read_to_string(filename).expect("Could not read file")
    };
    use std::path::Path;
    let output_dir = Path::new("tests/");
    if output_dir.exists() {
        std::fs::remove_dir_all(output_dir).ok();
    }
    std::fs::create_dir_all(output_dir).unwrap();
    runner::snp_diff_from_toml(&toml).unwrap();
    println!("Hello, world!");
}
