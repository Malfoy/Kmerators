use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

const FASTA_BUFFER_CAPACITY: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone, Default)]
pub struct RunReport {
    pub done: Vec<String>,
    pub failed: Vec<String>,
    pub multiple: Vec<String>,
    pub warning: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct InputMetadata {
    pub role: String,
    pub path: PathBuf,
    pub size_bytes: u64,
    pub modified_unix_seconds: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct PhaseTiming {
    pub label: String,
    pub elapsed: Duration,
}

#[derive(Debug, Clone)]
pub struct ReportMetadata {
    pub command: String,
    pub species: String,
    pub release: String,
    pub kmer_specs: Vec<(usize, usize)>,
    pub hash_table_count: usize,
    pub max_on_transcriptome: u32,
    pub max_on_genome: u32,
    pub threads: usize,
    pub stringent: bool,
    pub write_kmers: bool,
    pub inputs: Vec<InputMetadata>,
    pub phases: Vec<PhaseTiming>,
    pub total_elapsed: Duration,
    pub peak_rss_kib: Option<u64>,
}

pub struct LazyFastaWriter {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
}

impl LazyFastaWriter {
    pub fn new(path: PathBuf) -> Self {
        Self { path, writer: None }
    }

    pub fn write_record_with<F>(&mut self, write_header: F, seq: &[u8]) -> Result<()>
    where
        F: FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
    {
        let writer = self.ensure_writer()?;
        writer.write_all(b">")?;
        write_header(writer)?;
        writer.write_all(b"\n")?;
        writer.write_all(seq)?;
        writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        if let Some(writer) = &mut self.writer {
            writer
                .flush()
                .with_context(|| format!("failed to flush {}", self.path.display()))?;
        }
        Ok(())
    }

    fn ensure_writer(&mut self) -> Result<&mut BufWriter<File>> {
        if self.writer.is_none() {
            if let Some(parent) = self.path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }
            let file = File::create(&self.path)
                .with_context(|| format!("failed to create {}", self.path.display()))?;
            self.writer = Some(BufWriter::with_capacity(FASTA_BUFFER_CAPACITY, file));
        }
        Ok(self.writer.as_mut().expect("writer just initialized"))
    }
}

pub fn print_report(report: &RunReport) {
    if !report.done.is_empty() {
        println!("\nDone ({}):", report.done.len());
        for line in report.done.iter().take(15) {
            println!("  - {line}");
        }
        if report.done.len() > 15 {
            println!("  - ...");
        }
    }
    if !report.multiple.is_empty() {
        println!("\nMultiple responses ({}):", report.multiple.len());
        for line in report.multiple.iter().take(15) {
            println!("  - {line}");
        }
    }
    if !report.failed.is_empty() {
        println!("\nFailed ({}):", report.failed.len());
        for line in report.failed.iter().take(15) {
            println!("  - {line}");
        }
    }
    if !report.warning.is_empty() {
        println!("\nWarning ({}):", report.warning.len());
        for line in report.warning.iter().take(15) {
            println!("  - {line}");
        }
    }
}

pub fn write_markdown_report(
    path: &Path,
    metadata: &ReportMetadata,
    report: &RunReport,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let mut text = String::new();
    text.push_str("# kmerators report\n\n");
    text.push_str(&format!(
        "**Command:** `{}`\n\n",
        escape_inline_code(&metadata.command)
    ));
    text.push_str(&format!("**Version:** `{}`\n\n", env!("CARGO_PKG_VERSION")));
    text.push_str(&format!("**Species:** `{}`\n\n", metadata.species));
    text.push_str(&format!("**Release:** `{}`\n\n", metadata.release));
    text.push_str("## Parameters\n\n");
    text.push_str("| Setting | Value |\n|---|---:|\n");
    text.push_str(&format!(
        "| k/minimizer | {} |\n",
        metadata
            .kmer_specs
            .iter()
            .map(|(k, m)| format!("{k}/{m}"))
            .collect::<Vec<_>>()
            .join(", ")
    ));
    text.push_str(&format!(
        "| Minimizer hash tables | {} |\n",
        metadata.hash_table_count
    ));
    text.push_str(&format!(
        "| Maximum transcriptome count | {} |\n",
        metadata.max_on_transcriptome
    ));
    text.push_str(&format!(
        "| Maximum genome count | {} |\n",
        metadata.max_on_genome
    ));
    text.push_str(&format!("| Worker threads | {} |\n", metadata.threads));
    text.push_str(&format!(
        "| Stringent gene mode | {} |\n",
        metadata.stringent
    ));
    text.push_str(&format!("| Write k-mers | {} |\n\n", metadata.write_kmers));

    text.push_str("## Inputs\n\n");
    text.push_str("| Role | Path | Bytes | Modified (Unix seconds) |\n");
    text.push_str("|---|---|---:|---:|\n");
    for input in &metadata.inputs {
        let modified = input
            .modified_unix_seconds
            .map(|value| value.to_string())
            .unwrap_or_else(|| "unknown".to_string());
        text.push_str(&format!(
            "| {} | `{}` | {} | {} |\n",
            input.role,
            escape_inline_code(&input.path.display().to_string()),
            input.size_bytes,
            modified
        ));
    }
    text.push('\n');

    text.push_str("## Performance\n\n");
    text.push_str("| Phase | Elapsed seconds |\n|---|---:|\n");
    for phase in &metadata.phases {
        text.push_str(&format!(
            "| {} | {:.3} |\n",
            phase.label,
            phase.elapsed.as_secs_f64()
        ));
    }
    text.push_str(&format!(
        "| Total extraction | {:.3} |\n\n",
        metadata.total_elapsed.as_secs_f64()
    ));
    if let Some(peak_rss_kib) = metadata.peak_rss_kib {
        text.push_str(&format!("**Peak RSS:** `{peak_rss_kib} KiB`\n\n"));
    }

    push_section(&mut text, "Done", &report.done);
    push_section(&mut text, "Multiple responses", &report.multiple);
    push_section(&mut text, "Failed", &report.failed);
    push_section(&mut text, "Warning", &report.warning);
    std::fs::write(path, text).with_context(|| format!("failed to write {}", path.display()))
}

fn escape_inline_code(value: &str) -> String {
    value.replace('`', "\\`")
}

fn push_section(text: &mut String, title: &str, lines: &[String]) {
    text.push_str(&format!("## {title} ({})\n\n", lines.len()));
    if lines.is_empty() {
        text.push_str("None\n\n");
    } else {
        for line in lines {
            text.push_str(&format!("- {line}\n"));
        }
        text.push('\n');
    }
}
