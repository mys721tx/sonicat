use bio::io::fasta;
use clap::{Arg, Command, value_parser};
use rand::{SeedableRng, rngs::StdRng};
use sonicat::Mutator;
use std::fs::File;
use std::io::{self, BufWriter};
use std::process;

// from Brodin et al. 2013, doi:10.1371/journal.pone.0070388
const DEFAULT_SUBSTITUTION: f64 = 0.000057;
const DEFAULT_INSERTION: f64 = 0.000069;
const DEFAULT_DELETION: f64 = 0.0016;

fn main() {
    let matches = Command::new("Muta")
        .about("in silico mutation of FASTA sequences.")
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
            Arg::new("substitution")
                .short('s')
                .long("substitution")
                .value_name("SUBSTITUTION")
                .help(format!(
                    "Probability of substitution per nucleotide, default to {}",
                    DEFAULT_SUBSTITUTION
                ))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("insertion")
                .short('n')
                .long("insertion")
                .value_name("INSERTION")
                .help(format!(
                    "Probability of insertion per nucleotide, default to {}",
                    DEFAULT_INSERTION
                ))
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("deletion")
                .short('d')
                .long("deletion")
                .value_name("DELETION")
                .help(format!(
                    "Probability of deletion per nucleotide, default to {}",
                    DEFAULT_DELETION
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

    let substitution = matches
        .get_one::<f64>("substitution")
        .copied()
        .unwrap_or(DEFAULT_SUBSTITUTION);
    let insertion = matches
        .get_one::<f64>("insertion")
        .copied()
        .unwrap_or(DEFAULT_INSERTION);
    let deletion = matches
        .get_one::<f64>("deletion")
        .copied()
        .unwrap_or(DEFAULT_DELETION);

    let mutator = Mutator::new(substitution, insertion, deletion).unwrap_or_else(|e| {
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

    let mut buf = Vec::new();

    for record in reader.records() {
        let record = record.unwrap();

        buf.clear();
        mutator.mutate_into(record.seq(), &mut rng, &mut buf);
        writer.write(record.id(), record.desc(), &buf).unwrap();
    }
    writer.flush().unwrap();
}
