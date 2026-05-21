use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn tiny_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("tiny")
        .join(name)
}

fn read_text(path: PathBuf) -> String {
    std::fs::read_to_string(&path).unwrap_or_else(|err| {
        panic!("failed to read {}: {err}", path.display());
    })
}

fn run_tiny(out_dir: &Path, config_dir: &Path, extra_args: &[&str]) -> Output {
    run_tiny_with_paths(
        &tiny_fixture("query.fa"),
        &tiny_fixture("transcriptome.fa"),
        &tiny_fixture("genome.fa"),
        out_dir,
        config_dir,
        extra_args,
    )
}

fn run_tiny_with_paths(
    query: &Path,
    transcriptome: &Path,
    genome: &Path,
    out_dir: &Path,
    config_dir: &Path,
    extra_args: &[&str],
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_kmerators"));
    command
        .env("XDG_CONFIG_HOME", config_dir)
        .arg("-f")
        .arg(query)
        .arg("--transcriptome-fasta")
        .arg(transcriptome)
        .arg("-g")
        .arg(genome)
        .arg("-S")
        .arg("toy_species")
        .arg("-r")
        .arg("1")
        .arg("-k")
        .arg("5")
        .arg("-m")
        .arg("3")
        .arg("--hash-tables")
        .arg("64")
        .arg("-o")
        .arg(out_dir)
        .arg("-t")
        .arg("2")
        .arg("-y");
    command.args(extra_args);
    command.output().expect("failed to run kmerators")
}

fn write_gzip(path: &Path, bytes: &[u8]) {
    let file = std::fs::File::create(path).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap();
}

fn write_zstd(path: &Path, bytes: &[u8]) {
    let compressed = zstd::stream::encode_all(bytes, 0).unwrap();
    std::fs::write(path, compressed).unwrap();
}

fn write_xz(path: &Path, bytes: &[u8]) {
    let file = std::fs::File::create(path).unwrap();
    let mut encoder = liblzma::write::XzEncoder::new(file, 6);
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap();
}

fn assert_success(output: Output) {
    assert!(
        output.status.success(),
        "kmerators failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn tiny_fasta_cli_defaults_to_no_kmer_file() {
    let tmp = tempfile::tempdir().expect("failed to create tempdir");
    let out_dir = tmp.path().join("out");
    let config_dir = tmp.path().join("xdg");

    assert_success(run_tiny(&out_dir, &config_dir, &[]));

    assert!(!out_dir.join("kmers.fa").exists());
    for name in ["contigs.fa", "masked.fa"] {
        let expected = read_text(tiny_fixture(&format!("expected/{name}")));
        let actual = read_text(out_dir.join(name));
        assert_eq!(actual, expected, "{name} did not match the toy fixture");
    }

    let report = read_text(out_dir.join("report.md"));
    assert!(report.contains("q1: q1 - kmers/contigs: 2/2 (fasta)"));
}

#[test]
fn tiny_fasta_cli_writes_kmers_when_requested() {
    let tmp = tempfile::tempdir().expect("failed to create tempdir");
    let out_dir = tmp.path().join("out");
    let config_dir = tmp.path().join("xdg");

    assert_success(run_tiny(&out_dir, &config_dir, &["--write-kmers"]));

    for name in ["kmers.fa", "contigs.fa", "masked.fa"] {
        let expected = read_text(tiny_fixture(&format!("expected/{name}")));
        let actual = read_text(out_dir.join(name));
        assert_eq!(actual, expected, "{name} did not match the toy fixture");
    }
}

#[test]
fn tiny_fasta_cli_honors_transcriptome_threshold() {
    let tmp = tempfile::tempdir().expect("failed to create tempdir");
    let out_dir = tmp.path().join("out");
    let config_dir = tmp.path().join("xdg");

    assert_success(run_tiny(
        &out_dir,
        &config_dir,
        &["--write-kmers", "--max-on-transcriptome", "1"],
    ));

    assert_eq!(
        read_text(out_dir.join("kmers.fa")),
        concat!(
            ">q1.kmer2 ct:1\n",
            "CGTAC\n",
            ">q1.kmer3 ct:1\n",
            "GTACC\n",
            ">q1.kmer4 ct:1\n",
            "TACCC\n",
            ">q1.kmer5 ct:1\n",
            "ACCCC\n",
        )
    );
    assert_eq!(
        read_text(out_dir.join("contigs.fa")),
        ">q1.contig_1 (at position 2)\nCGTACCCC\n"
    );
    assert_eq!(
        read_text(out_dir.join("masked.fa")),
        ">q1.kmer1 genome:2 transcriptome:0\nACGTA\n"
    );

    let report = read_text(out_dir.join("report.md"));
    assert!(report.contains("q1: q1 - kmers/contigs: 4/1 (fasta)"));
}

#[test]
fn tiny_fasta_cli_repeated_k_writes_per_k_outputs() {
    let tmp = tempfile::tempdir().expect("failed to create tempdir");
    let out_dir = tmp.path().join("out");
    let config_dir = tmp.path().join("xdg");

    assert_success(run_tiny(
        &out_dir,
        &config_dir,
        &["--write-kmers", "-k", "6"],
    ));

    assert!(!out_dir.join("kmers.fa").exists());
    for name in ["kmers.fa", "contigs.fa", "masked.fa"] {
        let expected = read_text(tiny_fixture(&format!("expected/{name}")));
        let actual = read_text(out_dir.join("k5").join(name));
        assert_eq!(actual, expected, "k5/{name} did not match the toy fixture");
    }

    assert_eq!(
        read_text(out_dir.join("k6/kmers.fa")),
        concat!(
            ">q1.kmer1 ct:1\n",
            "ACGTAC\n",
            ">q1.kmer2 ct:1\n",
            "CGTACC\n",
            ">q1.kmer4 ct:2\n",
            "TACCCC\n",
        )
    );
    assert_eq!(
        read_text(out_dir.join("k6/contigs.fa")),
        concat!(
            ">q1.contig_1 (at position 1)\n",
            "ACGTACC\n",
            ">q1.contig_2 (at position 4)\n",
            "TACCCC\n",
        )
    );
    assert_eq!(
        read_text(out_dir.join("k6/masked.fa")),
        ">q1.kmer3 genome:0 transcriptome:1\nGTACCC\n"
    );

    let report = read_text(out_dir.join("report.md"));
    assert!(report.contains("k=5: q1: q1 - kmers/contigs: 2/2 (fasta)"));
    assert!(report.contains("k=6: q1: q1 - kmers/contigs: 3/2 (fasta)"));
}

#[test]
fn tiny_fasta_cli_accepts_gzip_zstd_and_xz_inputs() {
    let tmp = tempfile::tempdir().expect("failed to create tempdir");
    let out_dir = tmp.path().join("out");
    let config_dir = tmp.path().join("xdg");
    let query = tmp.path().join("query.fa.gz");
    let transcriptome = tmp.path().join("transcriptome.fa.zst");
    let genome = tmp.path().join("genome.fa.xz");

    write_gzip(
        &query,
        std::fs::read(tiny_fixture("query.fa")).unwrap().as_slice(),
    );
    write_zstd(
        &transcriptome,
        std::fs::read(tiny_fixture("transcriptome.fa"))
            .unwrap()
            .as_slice(),
    );
    write_xz(
        &genome,
        std::fs::read(tiny_fixture("genome.fa")).unwrap().as_slice(),
    );

    assert_success(run_tiny_with_paths(
        &query,
        &transcriptome,
        &genome,
        &out_dir,
        &config_dir,
        &["--write-kmers"],
    ));

    for name in ["kmers.fa", "contigs.fa", "masked.fa"] {
        let expected = read_text(tiny_fixture(&format!("expected/{name}")));
        let actual = read_text(out_dir.join(name));
        assert_eq!(actual, expected, "{name} did not match the toy fixture");
    }
}
