# kmerators-wasm

Browser-first WebAssembly prototype for local k-mer specificity filtering.

This directory is self-contained inside the KmeratoRS repository. It has its
own Rust workspace and web app so the browser prototype can evolve without
changing the native command-line crate layout.

## Scope

- Local FASTA/FASTQ query, transcriptome, and genome files.
- Transcriptome and genome references are optional filters.
- Optional default Ensembl human GRCh38 references can be streamed directly:
  genome primary assembly, transcriptome cDNA, and transcriptome ncRNA.
- Streaming reference scans so the genome is not loaded into memory.
- Plain, gzip, zstd, and xz input FASTA/FASTQ files.
- Plain, gzip, and zstd output downloads.
- Multiple k-mer sizes in one run.
- Exact k-mers up to 31 bases and hashed k-mer keys above 31 bases.
- Configurable minimizer length.
- Rust core compiled to `wasm32-unknown-unknown`.
- Browser UI that runs the core inside a Web Worker.

## Build

```sh
./scripts/build-wasm.sh
```

## Run

```sh
cd web
npm install
npm run start
```

Open the local URL printed by Vite.

## Test

```sh
cargo test
./scripts/build-wasm.sh
cd web && npm run check && npm run build
```

## Memory Notes

Reference FASTA/FASTQ files are decoded and scanned in chunks. They are not kept in
memory after scanning. Query sequences and query k-mer occurrences are kept in
memory because filtering and retained contig assembly need their original
coordinates and sequence text.

The Rust core still materializes final text outputs during finalization, but
the browser UI does not retain plain `contigs.fa`: the worker immediately
compresses it as zstd level `-1`, posts the compressed bytes to the UI, and
drops the plain text field. `kmers.fa` is never stored; plain and gzip downloads
are generated on demand from a streaming decompression of the stored contigs.
Zstd downloads for generated files use the current simple zstd compressor, so
that selected generated output is temporarily collected before compression.

The default remote references use Ensembl `current_fasta` HTTPS endpoints. They
require network access and are large, especially the human genome FASTA.
