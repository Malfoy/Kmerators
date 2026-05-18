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
READS="${READS:-250000}"
READ_LEN="${READ_LEN:-150}"
QUERY_LEN="${QUERY_LEN:-1000}"
K="${K:-31}"
THREADS="${THREADS:-$(nproc)}"
READ_MODE="${READ_MODE:-constant}"
M="${M:-9}"
HASH_TABLES="${HASH_TABLES:-1024}"

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
: > "${WORK}/timings.tsv"

time_cmd() {
  local label="$1"
  shift
  /usr/bin/time -f "${label}\t%e\t%M" -a -o "${WORK}/timings.tsv" "$@"
}

time_cmd data_generation "${VENV}/bin/python" - <<'PY' "${WORK}" "${READS}" "${READ_LEN}" "${QUERY_LEN}" "${READ_MODE}"
import pickle
import random
import sys
from pathlib import Path

work = Path(sys.argv[1])
reads = int(sys.argv[2])
read_len = int(sys.argv[3])
query_len = int(sys.argv[4])
read_mode = sys.argv[5]
rng = random.Random(20260513)
alphabet = "ACGT"

def randseq(n):
    return "".join(rng.choice(alphabet) for _ in range(n))

query = randseq(query_len)
transcriptome = "T" * (query_len + 200)

(work / "query.fa").write_text(f">q1\n{query}\n")
(work / "transcriptome.fa").write_text(f">TR1\n{transcriptome}\n")

def write_random_reads(out, start, count):
    for i in range(start, start + count):
        out.write(f">read_{i}\n{randseq(read_len)}\n")

def write_constant_reads(out, start, count):
    chunk_size = 1000
    seq = "T" * read_len
    written = 0
    while written < count:
        n = min(chunk_size, count - written)
        base = start + written
        out.write("".join(f">read_{base + j}\n{seq}\n" for j in range(n)))
        written += n

write_reads = write_constant_reads if read_mode == "constant" else write_random_reads

with (work / "reads.fa").open("w") as out:
    write_reads(out, 0, reads // 2)
    out.write(">query_once\n")
    out.write(query + "\n")
    write_reads(out, reads // 2, reads - reads // 2)

geneinfo = {
    "assembly": "ASM",
    "chr": ["reads"],
    "version": 1,
    "gene": {},
    "symbol": {},
    "alias": {},
    "transcript": {},
}
datadir = work / "py_data"
(datadir / "test_species.ASM.1.geneinfo.pkl").write_bytes(pickle.dumps(geneinfo))
(datadir / "test_species.ASM.1.transcriptome.pkl").write_bytes(pickle.dumps({"TR1": transcriptome}))
(datadir / "test_species.ASM.1.report.md").write_text("# large k31 fixture\n")
PY

echo "Dataset:"
echo "  reads: ${READS}"
echo "  read length: ${READ_LEN}"
echo "  read bases: $((READS * READ_LEN + QUERY_LEN))"
echo "  query length: ${QUERY_LEN}"
echo "  k: ${K}"
echo "  minimizer length: ${M}"
echo "  hash tables: ${HASH_TABLES}"
echo "  threads: ${THREADS}"
echo "  read mode: ${READ_MODE}"

time_cmd jellyfish_genome_index \
  "${JELLYFISH}" count --canonical -m "${K}" -s 100M -t "${THREADS}" \
  "${WORK}/reads.fa" -o "${WORK}/reads.k${K}.jf"

time_cmd jellyfish_transcriptome_index \
  "${JELLYFISH}" count -m "${K}" -s 100000 -t "${THREADS}" \
  "${WORK}/transcriptome.fa" -o "${WORK}/py_data/test_species.ASM.1.k${K}.transcriptome.jf"

time_cmd python_kmerator \
  env PYTHONPATH="${PY_KMERATOR_MODULES}" PATH="$(dirname "${JELLYFISH}"):${PATH}" \
  "${VENV}/bin/python" "${PY_KMERATOR}" \
  -f "${WORK}/query.fa" \
  -d "${WORK}/py_data" \
  -g "${WORK}/reads.k${K}.jf" \
  -S test_species \
  -r 1 \
  -k "${K}" \
  -o "${WORK}/py_out" \
  -t "${THREADS}" \
  -y \
  --keep > "${WORK}/python.stdout" 2> "${WORK}/python.stderr"

time_cmd rust_kmerator \
  "${RUST_KMERATORS}" \
  -f "${WORK}/query.fa" \
  --transcriptome-fasta "${WORK}/transcriptome.fa" \
  -g "${WORK}/reads.fa" \
  -S test_species \
  -r 1 \
  -k "${K}" \
  -m "${M}" \
  --hash-tables "${HASH_TABLES}" \
  -o "${WORK}/rust_out" \
  -t "${THREADS}" \
  -y > "${WORK}/rust.stdout" 2> "${WORK}/rust.stderr"

count_records() {
  local path="$1"
  grep -c '^>' "${path}"
}

sha() {
  sha256sum "$1" | awk '{print $1}'
}

py_kmers_sha="$(sha "${WORK}/py_out/kmers.fa")"
rs_kmers_sha="$(sha "${WORK}/rust_out/kmers.fa")"
py_contigs_sha="$(sha "${WORK}/py_out/contigs.fa")"
rs_contigs_sha="$(sha "${WORK}/rust_out/contigs.fa")"

echo
echo "Timings (label, wall_seconds, max_rss_kb):"
cat "${WORK}/timings.tsv"

echo
echo "Output summary:"
echo "  Python kmers:  $(count_records "${WORK}/py_out/kmers.fa") records, sha256=${py_kmers_sha}"
echo "  Rust kmers:    $(count_records "${WORK}/rust_out/kmers.fa") records, sha256=${rs_kmers_sha}"
echo "  Python contigs: $(count_records "${WORK}/py_out/contigs.fa") records, sha256=${py_contigs_sha}"
echo "  Rust contigs:   $(count_records "${WORK}/rust_out/contigs.fa") records, sha256=${rs_contigs_sha}"

if [[ "${py_kmers_sha}" != "${rs_kmers_sha}" ]]; then
  echo "kmers.fa differs" >&2
  diff -u "${WORK}/py_out/kmers.fa" "${WORK}/rust_out/kmers.fa" | head -200 >&2
  exit 1
fi

if [[ "${py_contigs_sha}" != "${rs_contigs_sha}" ]]; then
  echo "contigs.fa differs" >&2
  diff -u "${WORK}/py_out/contigs.fa" "${WORK}/rust_out/contigs.fa" | head -200 >&2
  exit 1
fi

echo
echo "Experiment passed. Work directory: ${WORK}"
