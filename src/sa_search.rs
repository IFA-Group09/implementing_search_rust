use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::str;
use clap::Parser;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};
use crate::fasta::Record;

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

fn naive_binary_search(reference: &str, read: &str, sa: &Vec<usize>) -> (i32, i32) {
    let mut min_index = 0;
    let mut max_index = sa.len();

    while min_index < max_index {
        let c = (min_index + max_index) / 2;

        if reference[sa[c]..] < *read {
            min_index = c + 1;
        } else {
            max_index = c;
        }
    }

    let first = min_index as i32;
    max_index = sa.len();
    while min_index < max_index {
        let c = (min_index + max_index) / 2;
        if *read < reference[sa[c]..] {
            max_index = c;
        } else {
            min_index = c + 1;
        }
    }

    let last = max_index as i32;
    if (first > last) || !(reference[sa[first as usize]..].starts_with(read)) {
        return (-1, -1);
    }
    (first, last)
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

    for reference in references.records() {
        let reference = str::from_utf8(reference.expect("Error reading reference sequence").seq())
            .unwrap()
            .to_string()
            + "$";
        let sa = suffix_array(reference.as_bytes());
        let start_time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs();

        let reads_reader = BufReader::new(GzDecoder::new(File::open(reads_path.clone()).unwrap()));
        let reads = fasta::Reader::new(reads_reader);
        let mut read_num = 0;
        for read in &reads_vec {
            let (first, last) =
                naive_binary_search(&reference, str::from_utf8(read.seq()).unwrap(), &sa);
            if first > 0 {
                for _ in 0..last - first + 1 {
                    println!("{}", str::from_utf8(read.seq()).unwrap());
                }
            }

            if read_num % 10 == 0 {
                writeln!(benchmark_file, "sa,{:?},{},{}", reads_path, SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs() - start_time, read_num);
            }
            read_num += 1;
        }
    }
}
