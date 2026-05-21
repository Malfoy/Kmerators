use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const FASTA_BUFFER_CAPACITY: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone, Default)]
pub struct RunReport {
    pub done: Vec<String>,
    pub failed: Vec<String>,
    pub multiple: Vec<String>,
    pub warning: Vec<String>,
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
    command: &str,
    species: &str,
    release: &str,
    report: &RunReport,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let mut text = String::new();
    text.push_str("# kmerators report\n\n");
    text.push_str(&format!("**Command:** `{command}`\n\n"));
    text.push_str(&format!("**Specie:** `{species}`\n\n"));
    text.push_str(&format!("**Release:** `{release}`\n\n"));
    push_section(&mut text, "Done", &report.done);
    push_section(&mut text, "Multiple responses", &report.multiple);
    push_section(&mut text, "Failed", &report.failed);
    push_section(&mut text, "Warning", &report.warning);
    std::fs::write(path, text).with_context(|| format!("failed to write {}", path.display()))
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
