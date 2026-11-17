use bio::io::fasta;
use clap::{Arg, Command, value_parser};
use rand::{
    Rng,
    distributions::{Uniform, WeightedIndex},
    rngs::ThreadRng,
    thread_rng,
};
use std::fs::File;
use std::io;

// from Brodin et al. 2013, doi:10.1371/journal.pone.0070388
const DEFAULT_SUBSTITUTION: f64 = 0.000057;
const DEFAULT_INSERTION: f64 = 0.000069;
const DEFAULT_DELETION: f64 = 0.0016;

struct Mutator {
    weights: [f64; 4],
    rng: ThreadRng,
}

impl Mutator {
    fn new(s: f64, i: f64, d: f64) -> Mutator {
        Mutator {
            weights: [s, i, d, 1.0 - s - i - d],
            rng: thread_rng(),
        }
    }
    fn mutate(&mut self, b: u8) -> Option<u8> {
        let dist = WeightedIndex::new(self.weights).unwrap();
        let fate = self.rng.sample(dist);
        match (fate, b) {
            (0, _) => Some({
                let dist = Uniform::from(0..3);
                let nuc = self.rng.sample(dist);
                match nuc {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => b,
                }
            }),
            (1, _) => None,
            (2, b'A') => Some({
                let dist = Uniform::from(0..2);
                let nuc = self.rng.sample(dist);
                match nuc {
                    0 => b'C',
                    1 => b'G',
                    2 => b'T',
                    _ => b,
                }
            }),
            (2, b'C') => Some({
                let dist = Uniform::from(0..2);
                let nuc = self.rng.sample(dist);
                match nuc {
                    0 => b'A',
                    1 => b'G',
                    2 => b'T',
                    _ => b,
                }
            }),
            (2, b'G') => Some({
                let dist = Uniform::from(0..2);
                let nuc = self.rng.sample(dist);
                match nuc {
                    0 => b'A',
                    1 => b'C',
                    2 => b'T',
                    _ => b,
                }
            }),
            (2, b'T') => Some({
                let dist = Uniform::from(0..2);
                let nuc = self.rng.sample(dist);
                match nuc {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b,
                }
            }),
            _ => Some(b),
        }
    }
}

fn main() {
    let matches = Command::new("matchesMuta")
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

    let mut mutator = Mutator::new(substitution, insertion, deletion);

    for record in reader.records() {
        let record = record.unwrap();

        let mut buf = Vec::with_capacity(record.seq().len() * 2);

        for r in record.seq().iter() {
            if let Some(x) = mutator.mutate(*r) {
                buf.push(x);
            }
        }
        writer.write(record.id(), record.desc(), &buf).unwrap();
    }
}
