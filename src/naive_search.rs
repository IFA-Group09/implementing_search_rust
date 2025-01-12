use bio::io::fasta;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::str;
use std::path::PathBuf;
use crate::fasta::Record;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about= None)]
struct Args {
    /// Path to reads in a gzipped fasta file 
    reads: PathBuf,
    /// Path to references in a gzipped fasta file
    references: PathBuf,
    /// Number of reads to use
    read_num: usize
}

fn naive_search(read: &str, reference: &str) {
    let matches: Vec<_> = reference.match_indices(&read).collect();
    for _m in matches {
        println!("{}", read);
    }
}

fn get_reads(reads_path: &PathBuf) -> Vec<Record> {
    let reads_reader = BufReader::new(GzDecoder::new(File::open(reads_path.clone()).unwrap()));
    let reads = fasta::Reader::new(reads_reader);
    let mut reads_vec: Vec<Record> = Vec::new();
    for read in reads.records() {
        let read = read.expect("error fetching read");
        reads_vec.push(read);
    }
    reads_vec
}

fn main() {
    let args = Args::parse();

    let reads_path = args.reads;
    let mut reads_vec: Vec<Record> = get_reads(&reads_path);

    while reads_vec.len() < args.read_num {
        let new_reads = get_reads(&reads_path);
        reads_vec.extend(new_reads);
    }
    reads_vec.truncate(args.read_num);

    let references_path = args.references;

    let references_reader = BufReader::new(GzDecoder::new(File::open(references_path).unwrap()));

    let references = fasta::Reader::new(references_reader);

    let mut benchmark_file = OpenOptions::new().create(true).append(true).open("rs_benchmark.csv").expect("error opening benchmark file");
    if benchmark_file.metadata().expect("error opening file metadata").len() == 0 {
        writeln!(benchmark_file, "method,reads_file,time,read_n");
    }

    let start_time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs();
    for reference in references.records() {
        let reference = reference.expect("error fetching reference");

        let reads_reader = BufReader::new(GzDecoder::new(File::open(reads_path.clone()).unwrap()));
        let reads = fasta::Reader::new(reads_reader);

        let mut read_num = 0;
        for read in &reads_vec {
            naive_search(
                str::from_utf8(read.seq()).unwrap(),
                str::from_utf8(reference.seq()).unwrap(),
            );

            if read_num % 10 == 0 {
                writeln!(benchmark_file, "naive,{:?},{},{}", reads_path, SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs() - start_time, read_num);
            }
            read_num += 1;
        }
    }
}
