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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::Path;

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

    #[test]
    fn normalize_dna_uppercases_ascii_bases() {
        assert_eq!(normalize_dna(b"acgtnN-"), b"ACGTNN-");
    }

    #[test]
    fn first_token_stops_at_ascii_whitespace() {
        assert_eq!(first_token(b"ENST0001.5 transcript name"), "ENST0001.5");
        assert_eq!(first_token(b"ENSG0001\tgene name"), "ENSG0001");
    }

    #[test]
    fn first_token_without_version_strips_suffix() {
        assert_eq!(
            first_token_without_version(b"ENST00000255409.4 description"),
            "ENST00000255409"
        );
    }

    #[test]
    fn write_fasta_record_uses_standard_two_line_format() {
        let mut out = Vec::new();

        write_fasta_record(&mut out, "query1", b"ACGT").unwrap();

        assert_eq!(out, b">query1\nACGT\n");
    }

    #[test]
    fn read_records_supports_gzip_zstd_and_xz_files() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = b">r1 description\nacgt\n>r2\nTTAA\n";

        let gz = tmp.path().join("records.fa.gz");
        let zst = tmp.path().join("records.fa.zst");
        let xz = tmp.path().join("records.fa.xz");
        write_gzip(&gz, fasta);
        write_zstd(&zst, fasta);
        write_xz(&xz, fasta);

        for path in [gz, zst, xz] {
            let records = read_records(&path).unwrap();
            assert_eq!(records.len(), 2, "failed for {}", path.display());
            assert_eq!(records[0].header, b"r1 description");
            assert_eq!(records[0].seq, b"ACGT");
            assert_eq!(records[1].header, b"r2");
            assert_eq!(records[1].seq, b"TTAA");
        }
    }
}
