use bio::io::fasta;
use clap::{App, Arg};
use rand::{thread_rng, Rng};
use rand_distr::Poisson;
use std::fs::File;
use std::io;

const DEFAULT_DEPTH: f64 = 50.0;
const DEFAULT_LENGTH: usize = 150;

fn main() {
    let matches = App::new("Sonicat")
        .about("in silico sonication of FASTA sequences.")
        .arg(
            Arg::with_name("in")
                .short("i")
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA file, default to stdin")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .value_name("OUTPUT")
                .help("Output FASTA file, default to stdout")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("depth")
                .short("d")
                .long("depth")
                .value_name("DEPTH")
                .help(format!("Average read depth, default to {}", DEFAULT_DEPTH).as_str())
                .takes_value(true),
        )
        .arg(
            Arg::with_name("length")
                .short("l")
                .long("length")
                .value_name("LENGTH")
                .help(format!("Average read length, default to {}", DEFAULT_LENGTH).as_str())
                .takes_value(true),
        )
        .get_matches();

    let fin: Box<dyn io::Read> = match matches.value_of("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };
    let reader = fasta::Reader::new(fin);

    let fout: Box<dyn io::Write> = match matches.value_of("out") {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut writer = fasta::Writer::new(fout);

    let depth = matches
        .value_of("depth")
        .map_or(DEFAULT_DEPTH, |x| x.parse().unwrap());
    let length = matches
        .value_of("length")
        .map_or(DEFAULT_LENGTH, |x| x.parse().unwrap());

    let poi = Poisson::new(depth).unwrap();

    let mut count = 0;

    for record in reader.records() {
        let record = record.unwrap();
        let record = record.seq().windows(length);

        for r in record {
            let v = thread_rng().sample(poi) as u64;
            for _ in 0..v {
                count += 1;

                let name = format!("seq_{}", count);

                writer.write(name.as_str(), Option::None, r).unwrap();
            }
        }
    }
}
