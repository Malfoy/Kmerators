# Large k=31 Python/Rust Experiment

This experiment compares Python `kmerator` and the Rust `kmerators` clone on a larger synthetic `--fasta-file` case:

- `k = 31`
- query length: 1000 bp
- default read dataset: 250,000 reads x 150 bp = 37.5 Mbp, plus one 1000 bp query-containing read
- transcriptome: unrelated local sequence, so retained k-mers are controlled by genome counts
- Rust minimizer routing: `m=9`, `1024` hash tables by default

Run:

```sh
./experiments/large-k31/run.sh
```

You can scale the read dataset:

```sh
READS=1000000 READ_LEN=150 ./experiments/large-k31/run.sh
```

You can also tune minimizer routing:

```sh
M=9 HASH_TABLES=1024 ./experiments/large-k31/run.sh
```

The script reports timings for:

- Jellyfish genome index build
- Jellyfish transcriptome index build
- Python `kmerator` extraction
- Rust `kmerators` extraction

It then compares `kmers.fa` and `contigs.fa` by SHA-256 and record counts.
