#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ATRAVERS_DIR="$(cd "${ROOT_DIR}/.." && pwd)"

PY_KMERATOR="${ATRAVERS_DIR}/kmerator/kmerator/kmerator.py"
PY_KMERATOR_MODULES="${ATRAVERS_DIR}/kmerator/kmerator"
RUST_KMERATORS="${ROOT_DIR}/target/release/kmerators"
JELLYFISH="${JELLYFISH:-/home/nadine/Code/Jellyfish/bin/jellyfish}"
VENV="${ROOT_DIR}/.venv-parity"
WORK="${SCRIPT_DIR}/work"

if [[ ! -x "${JELLYFISH}" ]]; then
  echo "Jellyfish binary not found: ${JELLYFISH}" >&2
  exit 1
fi

if [[ ! -x "${RUST_KMERATORS}" ]]; then
  cargo build --release --bin kmerators --manifest-path "${ROOT_DIR}/Cargo.toml"
fi

if [[ ! -x "${VENV}/bin/python" ]]; then
  python3 -m venv "${VENV}"
  "${VENV}/bin/python" -m pip install --upgrade pip
  "${VENV}/bin/python" -m pip install bs4 lxml requests
fi

rm -rf "${WORK}"
mkdir -p "${WORK}/py_data" "${WORK}/py_out" "${WORK}/rust_out"

cat > "${WORK}/genome.fa" <<'FA'
>chr1
TTTACGTAGGGACCCC
FA

cat > "${WORK}/transcriptome.fa" <<'FA'
>TR1
GGGGTACCCAAAA
FA

cat > "${WORK}/query.fa" <<'FA'
>q1
ACGTACCCC
FA

"${JELLYFISH}" count --canonical -m 5 -s 100000 -t 1 \
  "${WORK}/genome.fa" -o "${WORK}/genome.k5.jf"
"${JELLYFISH}" count -m 5 -s 100000 -t 1 \
  "${WORK}/transcriptome.fa" -o "${WORK}/py_data/test_species.ASM.1.k5.transcriptome.jf"

"${VENV}/bin/python" - <<'PY' "${WORK}/py_data"
import pickle
import sys
from pathlib import Path

datadir = Path(sys.argv[1])
geneinfo = {
    "assembly": "ASM",
    "chr": ["chr1"],
    "version": 1,
    "gene": {},
    "symbol": {},
    "alias": {},
    "transcript": {},
}
transcriptome = {"TR1": "GGGGTACCCAAAA"}
(datadir / "test_species.ASM.1.geneinfo.pkl").write_bytes(pickle.dumps(geneinfo))
(datadir / "test_species.ASM.1.transcriptome.pkl").write_bytes(pickle.dumps(transcriptome))
(datadir / "test_species.ASM.1.report.md").write_text("# tiny fixture\n")
PY

PYTHONPATH="${PY_KMERATOR_MODULES}" PATH="$(dirname "${JELLYFISH}"):${PATH}" \
  "${VENV}/bin/python" "${PY_KMERATOR}" \
  -f "${WORK}/query.fa" \
  -d "${WORK}/py_data" \
  -g "${WORK}/genome.k5.jf" \
  -S test_species \
  -r 1 \
  -k 5 \
  -o "${WORK}/py_out" \
  -t 1 \
  -y \
  --keep > "${WORK}/python.stdout" 2> "${WORK}/python.stderr"

"${RUST_KMERATORS}" \
  -f "${WORK}/query.fa" \
  --transcriptome-fasta "${WORK}/transcriptome.fa" \
  -g "${WORK}/genome.fa" \
  -S test_species \
  -r 1 \
  -k 5 \
  --hash-tables 1024 \
  -o "${WORK}/rust_out" \
  -t 2 \
  -y > "${WORK}/rust.stdout" 2> "${WORK}/rust.stderr"

echo "Python kmerator kmers:"
cat "${WORK}/py_out/kmers.fa"
echo
echo "Rust kmerators kmers:"
cat "${WORK}/rust_out/kmers.fa"
echo

echo "Diff kmers.fa"
diff -u "${WORK}/py_out/kmers.fa" "${WORK}/rust_out/kmers.fa"

echo "Diff contigs.fa"
diff -u "${WORK}/py_out/contigs.fa" "${WORK}/rust_out/contigs.fa"

echo "Masked outputs are not expected to be byte-identical; Python uses contig ids in masked FASTA headers for --fasta-file."
echo "Experiment passed. Work directory: ${WORK}"
