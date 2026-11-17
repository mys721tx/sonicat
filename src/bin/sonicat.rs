use bio::io::fasta;
use clap::{Arg, Command, value_parser};
use rand::{SeedableRng, rngs::StdRng};
use sonicat::Sonicator;
use std::fs::File;
use std::io::{self, BufWriter};
use std::process;

const DEFAULT_DEPTH: f64 = 50.0;
const DEFAULT_LENGTH: f64 = 150.0;
const DEFAULT_SD: f64 = 15.0;

fn main() {
    let matches = Command::new("Sonicat")
        .about("in silico sonication of FASTA sequences.")
        .after_help(
            "Each fragment is named seq_N and described as ID:START-END, \
             the 1-based inclusive region of source record ID it covers.",
        )
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
                .help(format!(
                    "Average per-base coverage, default to {}",
                    DEFAULT_DEPTH
                ))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("length")
                .short('l')
                .long("length")
                .value_name("LENGTH")
                .help(format!(
                    "Average fragment length, default to {}",
                    DEFAULT_LENGTH
                ))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("sd")
                .short('s')
                .long("sd")
                .value_name("SD")
                .help(format!(
                    "Standard deviation of fragment length, default to {}",
                    DEFAULT_SD
                ))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .value_name("SEED")
                .help("Seed for the random number generator, default to a random seed")
                .value_parser(value_parser!(u64)),
        )
        .get_matches();

    let depth = matches
        .get_one::<f64>("depth")
        .copied()
        .unwrap_or(DEFAULT_DEPTH);
    let length = matches
        .get_one::<f64>("length")
        .copied()
        .unwrap_or(DEFAULT_LENGTH);
    let sd = matches.get_one::<f64>("sd").copied().unwrap_or(DEFAULT_SD);

    let sonicator = Sonicator::new(depth, length, sd).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        process::exit(2);
    });

    let mut rng = match matches.get_one::<u64>("seed") {
        Some(&seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_rng(&mut rand::rng()),
    };

    let fin: Box<dyn io::Read> = match matches.get_one::<String>("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };
    let reader = fasta::Reader::new(fin);

    let fout: Box<dyn io::Write> = match matches.get_one::<String>("out") {
        Some(f) => Box::new(BufWriter::new(File::create(f).unwrap())),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut writer = fasta::Writer::new(fout);

    let mut count = 0;

    for record in reader.records() {
        let record = record.unwrap();
        let seq = record.seq();

        for fragment in sonicator.fragments(seq.len(), &mut rng) {
            count += 1;

            let name = format!("seq_{}", count);
            let desc = format!("{}:{}-{}", record.id(), fragment.start + 1, fragment.end);

            writer
                .write(&name, Some(&desc), fragment.slice(seq))
                .unwrap();
        }
    }
    writer.flush().unwrap();
}
