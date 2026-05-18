# Tiny Offline Example

This directory contains a minimal FASTA-only dataset for testing `kmerator-rs`
without Ensembl, Jellyfish, or the Python implementation.

Input:

- `query.fa`: one query sequence, `ACGTACCCC`
- `genome.fa`: one local reference sequence, `TTTACGTAGGGACCCC`
- `transcriptome.fa`: one local transcriptome sequence, `GGGGTACCCAAAA`

Run from `kmerator-rs/`:

```sh
bash examples/tiny/run.sh
```

The script writes to `examples/tiny/out` and compares the generated FASTA files
with `examples/tiny/expected`. Use `OUT_DIR=/tmp/kmerator-tiny-out bash
examples/tiny/run.sh` to keep generated files outside the repository.

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

The `-r 1` value is intentional: it keeps this local example fully offline by
avoiding the default `last` Ensembl release lookup.

Expected retained k-mers:

- `CGTAC`
- `ACCCC`

Expected contigs:

- `CGTAC`
- `ACCCC`

The report should contain `q1: q1 - kmers/contigs: 2/2 (fasta)`.

Multi-k smoke command:

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

This writes per-k outputs in `examples/tiny/out-multi/k5/` and
`examples/tiny/out-multi/k6/`, plus a combined report at
`examples/tiny/out-multi/report.md`.
