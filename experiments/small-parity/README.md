# Small Python/Rust Parity Experiment

This fixture compares Python `kmerator` and the Rust `kmerators` clone on one tiny `--fasta-file` query.

The synthetic data uses `k=5`:

- query: `ACGTACCCC`
- genome: `TTTACGTAGGGACCCC`
- transcriptome: `GGGGTACCCAAAA`

Expected retained k-mers:

- `CGTAC`
- `ACCCC`

Expected contigs:

- `CGTAC`
- `ACCCC`

Run:

```sh
./experiments/small-parity/run.sh
```

The script creates a fake Python kmerator dataset, builds Jellyfish indexes for Python, runs the Rust `kmerators` clone with `--transcriptome-fasta`, and diffs `kmers.fa` and `contigs.fa`.
