#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/out}"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

cargo run --release --bin kmerators --manifest-path "${ROOT_DIR}/Cargo.toml" -- \
  -f "${SCRIPT_DIR}/query.fa" \
  --transcriptome-fasta "${SCRIPT_DIR}/transcriptome.fa" \
  -g "${SCRIPT_DIR}/genome.fa" \
  -S toy_species \
  -r 1 \
  -k 5 \
  -m 3 \
  --hash-tables 64 \
  --write-kmers \
  -o "${OUT_DIR}" \
  -t 2 \
  -y

diff -u "${SCRIPT_DIR}/expected/kmers.fa" "${OUT_DIR}/kmers.fa"
diff -u "${SCRIPT_DIR}/expected/contigs.fa" "${OUT_DIR}/contigs.fa"
diff -u "${SCRIPT_DIR}/expected/masked.fa" "${OUT_DIR}/masked.fa"
grep -F "q1: q1 - kmers/contigs: 2/2 (fasta)" "${OUT_DIR}/report.md" >/dev/null

echo "Tiny example passed. Output directory: ${OUT_DIR}"
