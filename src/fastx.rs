use anyhow::{Context, Result};
use helicase::input::*;
use helicase::*;
use std::path::Path;

const FASTX_CONFIG: Config = ParserOptions::default().config();

#[derive(Debug, Clone)]
pub struct Record {
    pub header: Vec<u8>,
    pub seq: Vec<u8>,
}

pub fn read_records(path: &Path) -> Result<Vec<Record>> {
    let mut parser = FastxParser::<FASTX_CONFIG>::from_file(path)
        .with_context(|| format!("failed to open FASTA/FASTQ {}", path.display()))?;
    let mut records = Vec::new();
    while let Some(_) = parser.next() {
        records.push(Record {
            header: parser.get_header().to_vec(),
            seq: normalize_dna(parser.get_dna_string()),
        });
    }
    Ok(records)
}

pub fn normalize_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'a'..=b'z' => b - 32,
            _ => b,
        })
        .collect()
}

pub fn first_token(header: &[u8]) -> String {
    let token = header
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or(header);
    String::from_utf8_lossy(token).into_owned()
}

pub fn first_token_without_version(header: &[u8]) -> String {
    let token = first_token(header);
    token.split('.').next().unwrap_or(&token).to_string()
}

pub fn write_fasta_record<W: std::io::Write>(
    mut writer: W,
    header: &str,
    seq: &[u8],
) -> Result<()> {
    writer.write_all(b">")?;
    writer.write_all(header.as_bytes())?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}
