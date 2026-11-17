use bio::io::fasta;
use clap::{value_parser, Arg, Command};
use rand::{thread_rng, Rng};
use rand_distr::Poisson;
use std::fs::File;
use std::io;

const DEFAULT_DEPTH: f64 = 50.0;
const DEFAULT_LENGTH: usize = 150;

fn main() {
    let matches = Command::new("Sonicat")
        .about("in silico sonication of FASTA sequences.")
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output FASTA file, default to stdout"),
        )
        .arg(
            Arg::new("depth")
                .short('d')
                .long("depth")
                .value_name("DEPTH")
                .help(format!("Average read depth, default to {}", DEFAULT_DEPTH))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("length")
                .short('l')
                .long("length")
                .value_name("LENGTH")
                .help(format!(
                    "Average read length, default to {}",
                    DEFAULT_LENGTH
                ))
                .value_parser(value_parser!(usize)),
        )
        .get_matches();

    let fin: Box<dyn io::Read> = match matches.get_one::<String>("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };
    let reader = fasta::Reader::new(fin);

    let fout: Box<dyn io::Write> = match matches.get_one::<String>("out") {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut writer = fasta::Writer::new(fout);

    let depth = matches
        .get_one::<f64>("depth")
        .copied()
        .unwrap_or(DEFAULT_DEPTH);
    let length = matches
        .get_one::<usize>("length")
        .copied()
        .unwrap_or(DEFAULT_LENGTH);

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
