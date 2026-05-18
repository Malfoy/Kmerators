# kmerator-rs

`kmerator-rs` is a Rust implementation of `kmerator` for finding k-mers from
genes, transcripts, or local FASTA/FASTQ queries that are specific under genome
and transcriptome count thresholds.

The Rust version keeps the Python CLI shape and output style, but it does not
require Jellyfish for extraction. Genome input is read directly from
FASTA/FASTQ, optionally compressed. Query k-mers use exact `u64` keys for
`k <= 31`; larger k values use fixed `u64` hashes. K-mers are routed by
minimizer into hash tables for parallel counting.

## Requirements

- Rust toolchain with Cargo.
- A reference genome FASTA/FASTQ for extraction runs; gzip-compressed input is
  accepted.
- Internet access only for Ensembl-backed commands that use `-r last`,
  `--mk-dataset`, `--update-dataset`, `--last-avail`, or `--info`.

Jellyfish is not required for the Rust extraction path.

## Quick Start With Toy Data

The repository includes a tiny fully offline dataset in `examples/tiny/`.
It is the fastest way to check that the tool builds, runs, and produces stable
outputs on your machine:

```sh
bash examples/tiny/run.sh
```

The script runs the command below and compares `kmers.fa`, `contigs.fa`, and
`masked.fa` against checked-in expected files. It writes to
`examples/tiny/out` by default; set `OUT_DIR=/tmp/kmerator-tiny-out` to use a
throwaway output directory.

Equivalent manual command:

```sh
cargo run --release --bin kmerators -- \
  -f examples/tiny/query.fa \
  --transcriptome-fasta examples/tiny/transcriptome.fa \
  -g examples/tiny/genome.fa \
  -S toy_species \
  -r 1 \
  -k 5 \
  -m 3 \
  --hash-tables 64 \
  -o examples/tiny/out \
  -t 2 \
  -y
```

Expected output:

- `examples/tiny/out/kmers.fa` contains `CGTAC` and `ACCCC`
- `examples/tiny/out/contigs.fa` contains two contigs, `CGTAC` and `ACCCC`
- `examples/tiny/out/masked.fa` contains the rejected query k-mers
- `examples/tiny/out/report.md` summarizes the run

The `-r 1` argument keeps this example offline. If `-r` is omitted, the default
release is `last`, which asks Ensembl for the current release.

The same extraction can be run for several k-mer sizes in one command:

```sh
cargo run --release --bin kmerators -- \
  -f examples/tiny/query.fa \
  --transcriptome-fasta examples/tiny/transcriptome.fa \
  -g examples/tiny/genome.fa \
  -S toy_species \
  -r 1 \
  -k 5 \
  -k 6 \
  -m 3 \
  --hash-tables 64 \
  -o examples/tiny/out-multi \
  -t 2 \
  -y
```

For repeated `-k`, kmerator-rs loads the dataset and query sequences once, then
streams the transcriptome and genome once per counting phase while updating all
requested k-mer indexes. Outputs are written under one subdirectory per k value,
for example `out-multi/k5/` and `out-multi/k6/`, with a combined
`out-multi/report.md`.

## Build

Portable optimized build:

```sh
cargo build --release --bin kmerators
```

Host-optimized build:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release --no-default-features --bin kmerators
```

The built binary is `target/release/kmerators`.

## Input Modes

### Local FASTA/FASTQ Queries

Use `--fasta-file` for unannotated sequences. A local transcriptome can be
provided with `--transcriptome-fasta`, which bypasses kmerator dataset lookup:

```sh
target/release/kmerators \
  -f queries.fa \
  --transcriptome-fasta transcriptome.fa.gz \
  -g genome.fa.gz \
  -S my_species \
  -r 1 \
  -k 31 \
  -o output \
  -t 16 \
  -y
```

For FASTA/FASTQ queries, a k-mer is retained when:

- its transcriptome count is `<= --max-on-transcriptome` (`0` by default)
- its genome count is `<= --max-on-genome` (`1` by default)

### Ensembl Genes Or Transcripts

Build an Ensembl-backed dataset:

```sh
target/release/kmerators --mk-dataset -d data -S human -r last -y
```

Find specific k-mers for gene symbols, Ensembl gene IDs, aliases, or transcript
IDs:

```sh
target/release/kmerators \
  -s NPM1 ENST00000255409 \
  -d data \
  -g GRCh38.fa.gz \
  -o output \
  -t 16 \
  -y
```

You can also pass one selection file to `-s`; whitespace-separated tokens are
read and text after `#` is treated as a comment.

For gene mode, `--stringent` keeps only k-mers found in all isoforms. Without
`--stringent`, gene k-mers are retained when their transcriptome and isoform
presence counts are consistent and their genome count is at most one.

`--transcriptome-fasta` currently supports `--fasta-file` extraction only, not
gene/transcript selection.

## Outputs

Each extraction run writes:

- `kmers.fa`: retained specific k-mers
- `contigs.fa`: adjacent retained k-mers merged into contigs
- `masked.fa`: rejected query k-mers with genome/transcriptome counts
- `report.md`: done, failed, multiple-match, and warning summaries

## Useful Options

```sh
target/release/kmerators ... \
  -k 31 \
  -k 41 \
  -m 9 \
  --hash-tables 1024 \
  --max-on-transcriptome 0 \
  --max-on-genome 1 \
  -t 16
```

- `-k, --kmer-length`: k-mer length; repeatable; default `31`
- `-m, --minimizer-length`: minimizer length; default `min(9, k)` for each k
- `--hash-tables`: minimizer-routed hash table count; default `1024`
- `-t, --thread`: worker thread count; default is available CPU count
- `--tmpdir`: temporary directory
- `--keep`: keep intermediate files where applicable
- `-D, --debug`: print extra run details

## Tests And Fixtures

The checked-in tests are offline and deterministic. They cover low-level k-mer
helpers, local FASTA extraction, CLI validation for the offline transcriptome
path, and the tiny CLI golden fixture.

Run them with:

```sh
cargo test
```

Run the standalone toy smoke test:

```sh
bash examples/tiny/run.sh
```

Additional Python/Rust parity experiments live under `experiments/`. They are
useful for compatibility checks, but they require the Python `kmerator`
implementation and a Jellyfish binary.
