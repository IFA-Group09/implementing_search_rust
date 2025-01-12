Naive search can be run with:
```
cargo run --bin naive_search REFERENCES READS READ_N
```

where `REFERENCES` is a path to a gzipped FASTA containing references, `READS`  is a pth to a gzipped FASTA file containing reads, and `READ_N` is the number of reads to use.

Suffix array search can be run with:
```
cargo run --bin sa_search REFERENCES READS READ_N
```

where `REFERENCES` is a path to a gzipped FASTA containing references, `READS`  is a pth to a gzipped FASTA file containing reads, and `READ_N` is the number of reads to use.

