use anyhow::{Context, Result};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct RunReport {
    pub done: Vec<String>,
    pub failed: Vec<String>,
    pub multiple: Vec<String>,
    pub warning: Vec<String>,
}

pub fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for (header, seq) in records {
        writer.write_all(b">")?;
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(seq)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
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
