# kmerator-rs

Find specific k-mers from genes, transcripts, or local FASTA/FASTQ queries
without Jellyfish.

## Install

You need Rust with Cargo.

Install directly from GitHub:

```sh
cargo install --git https://github.com/Malfoy/Kmerators.git --locked --bin kmerators
```

Check the install:

```sh
kmerators --help
```

If `kmerators` is not found, add Cargo's binary directory to your `PATH`:

```sh
export PATH="$HOME/.cargo/bin:$PATH"
```

To run the examples or work on the code, clone the repository:

```sh
git clone https://github.com/Malfoy/Kmerators.git
cd Kmerators
cargo install --path . --locked --bin kmerators
```

## Run The Tiny Example

The tiny example is included in the repository, so run it from a clone:

```sh
cd Kmerators
bash examples/tiny/run.sh
```

If you installed directly with `cargo install --git`, clone the repository first
to get the example files.

The command should end with:

```text
Tiny example passed.
```

The tiny example is fully offline. It builds the tool, runs it on checked-in
FASTA files, and verifies `kmers.fa`, `contigs.fa`, and `masked.fa` against
expected outputs.

## What Files Do I Need?

For local FASTA/FASTQ mode, provide:

- `queries.fa`: sequences to extract specific k-mers from
- `transcriptome.fa`: transcriptome sequences used to reject non-specific k-mers
- `genome.fa`: reference genome sequences used to reject repeated k-mers

Any of these can be plain text or compressed as `.gz`, `.zst`, or `.xz`.

## Use Your Own FASTA/FASTQ

```sh
kmerators \
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

The `-r 1` value keeps this local-file command offline. If `-r` is omitted, the
default release is `last`, which asks Ensembl for the current release.

FASTA/FASTQ inputs can be plain text or compressed as `.gz`, `.zst`, or `.xz`.
Parsing is handled by `helicase` in all cases.

Results are written to `output/`:

- `kmers.fa`: retained specific k-mers
- `contigs.fa`: adjacent retained k-mers merged into contigs
- `masked.fa`: rejected query k-mers with genome/transcriptome counts
- `report.md`: run summary

Start by checking `output/report.md`; it lists completed, failed, ambiguous, and
warning cases for the run.

For FASTA/FASTQ queries, a k-mer is retained when its transcriptome count is
`<= --max-on-transcriptome` and its genome count is `<= --max-on-genome`.
Defaults are `0` and `1`.

## Use Several K-mer Sizes

Repeat `-k`:

```sh
kmerators \
  -f queries.fa \
  --transcriptome-fasta transcriptome.fa.gz \
  -g genome.fa.gz \
  -S my_species \
  -r 1 \
  -k 31 \
  -k 41 \
  -k 51 \
  -o output \
  -t 16 \
  -y
```

Outputs are written under one directory per size, for example `output/k31/`,
`output/k41/`, and `output/k51/`, with a combined `output/report.md`.

## Use Ensembl Genes Or Transcripts

Build a local Ensembl-backed dataset:

```sh
kmerators --mk-dataset -d data -S human -r last -y
```

This uses the network when `-r last` is used, and dataset creation can take time
because transcriptome and gene metadata are downloaded from Ensembl.

Find specific k-mers for gene symbols, Ensembl gene IDs, aliases, or transcript
IDs:

```sh
kmerators \
  -s NPM1 ENST00000255409 \
  -d data \
  -g GRCh38.fa.gz \
  -o output \
  -t 16 \
  -y
```

You can also pass one selection file to `-s`; whitespace-separated tokens are
read and text after `#` is treated as a comment.

## Requirements

- Rust toolchain with Cargo.
- A reference genome FASTA/FASTQ for extraction runs. Plain, `.gz`, `.zst`, and
  `.xz` inputs are accepted.
- Internet access only for Ensembl-backed commands that use `-r last`,
  `--mk-dataset`, `--update-dataset`, `--last-avail`, or `--info`.

Jellyfish is not required for the Rust extraction path.

## Build Without Installing

Portable optimized build:

```sh
cargo build --release --bin kmerators
```

The built binary is:

```sh
target/release/kmerators
```

Host-optimized build:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release --no-default-features --bin kmerators
```

## Useful Options

```sh
kmerators ... \
  -k 31 \
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

For gene mode, `--stringent` keeps only k-mers found in all isoforms. Without
`--stringent`, gene k-mers are retained when their transcriptome and isoform
presence counts are consistent and their genome count is at most one.

`--transcriptome-fasta` currently supports `--fasta-file` extraction only, not
gene/transcript selection.

## Tests

Run the full offline test suite:

```sh
cargo test --locked
```

Run the standalone toy smoke test:

```sh
bash examples/tiny/run.sh
```

Additional Python/Rust parity experiments live under `experiments/`. They are
useful for compatibility checks, but they require the Python `kmerator`
implementation and a Jellyfish binary.

## Implementation Notes

The Rust version keeps the Python CLI shape and output style. FASTA/FASTQ input
is parsed with `helicase`; plain, gzip, zstd, and xz streams are detected
automatically. Query k-mers use exact `u64` keys for `k <= 31`; larger k values
use fixed `u64` hashes. K-mers are routed by minimizer into hash tables for
parallel counting.
