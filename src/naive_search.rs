use bio::io::fasta;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::str;
use std::path::PathBuf;

use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about= None)]
struct Args {
    /// Path to reads in a gzipped fasta file 
    reads: PathBuf,
    /// Path to references in a gzipped fasta file
    references: PathBuf
}

fn naive_search(read: &str, reference: &str) {
    let matches: Vec<_> = reference.match_indices(&read).collect();
    for _m in matches {
        println!("{}", read);
    }
}

fn main() {
    let args = Args::parse();

    let reads_path = args.reads;
    let references_path = args.references;

    let references_reader = BufReader::new(GzDecoder::new(File::open(references_path).unwrap()));

    let references = fasta::Reader::new(references_reader);

    for reference in references.records() {
        let reference = reference.expect("error fetching reference");

        let reads_reader = BufReader::new(GzDecoder::new(File::open(reads_path.clone()).unwrap()));
        let reads = fasta::Reader::new(reads_reader);
        for read in reads.records() {
            let read = read.expect("error fetching read");
            naive_search(
                str::from_utf8(read.seq()).unwrap(),
                str::from_utf8(reference.seq()).unwrap(),
            );
        }
    }
}
