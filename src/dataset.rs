use crate::cli::Args;
use crate::ensembl;
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub const DATASET_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Dataset {
    pub species: String,
    pub release: String,
    pub assembly: String,
    pub geneinfo: GeneInfo,
    pub transcriptome: Transcriptome,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneInfo {
    pub assembly: String,
    pub chromosomes: Vec<String>,
    pub version: u32,
    pub genes: HashMap<String, Gene>,
    pub symbols: HashMap<String, Vec<String>>,
    pub aliases: HashMap<String, Vec<String>>,
    pub transcripts: HashMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub symbol: Option<String>,
    pub aliases: Vec<String>,
    pub canonical: Option<String>,
    pub transcripts: Vec<String>,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    pub desc: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcriptome {
    pub records: Vec<TranscriptRecord>,
    pub by_id: HashMap<String, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

impl Transcriptome {
    pub fn new(records: Vec<TranscriptRecord>) -> Self {
        let by_id = records
            .iter()
            .enumerate()
            .map(|(idx, rec)| (rec.id.to_ascii_uppercase(), idx))
            .collect();
        Self { records, by_id }
    }

    pub fn get(&self, id: &str) -> Option<&[u8]> {
        self.by_id
            .get(&id.to_ascii_uppercase())
            .and_then(|&idx| self.records.get(idx))
            .map(|rec| rec.seq.as_slice())
    }
}

pub fn dataset_basename(species: &str, assembly: &str, release: &str) -> String {
    format!("{species}.{assembly}.{release}")
}

pub fn dataset_path(args: &Args) -> Result<PathBuf> {
    let assembly = find_assembly_for_release(&args.datadir, &args.specie, &args.release)?
        .with_context(|| {
            format!(
                "dataset not found for species {} release {} in {}",
                args.specie,
                args.release,
                args.datadir.display()
            )
        })?;
    Ok(args.datadir.join(format!(
        "{}.dataset.bin",
        dataset_basename(&args.specie, &assembly, &args.release)
    )))
}

pub fn load_dataset(args: &Args) -> Result<Dataset> {
    let path = dataset_path(args)?;
    let bytes =
        std::fs::read(&path).with_context(|| format!("failed to read {}", path.display()))?;
    let (dataset, _): (Dataset, usize) =
        bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
            .with_context(|| format!("failed to decode {}", path.display()))?;
    if dataset.geneinfo.version != DATASET_VERSION {
        bail!(
            "dataset version {} is not supported by this kmerator-rs build",
            dataset.geneinfo.version
        );
    }
    Ok(dataset)
}

pub fn save_dataset(args: &Args, dataset: &Dataset) -> Result<PathBuf> {
    std::fs::create_dir_all(&args.datadir)
        .with_context(|| format!("failed to create {}", args.datadir.display()))?;
    let base = dataset_basename(&dataset.species, &dataset.assembly, &dataset.release);
    let path = args.datadir.join(format!("{base}.dataset.bin"));
    let bytes = bincode::serde::encode_to_vec(dataset, bincode::config::standard())?;
    std::fs::write(&path, bytes).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(path)
}

pub fn make_dataset(args: &Args) -> Result<()> {
    if !args.datadir.exists() {
        std::fs::create_dir_all(&args.datadir)
            .with_context(|| format!("failed to create {}", args.datadir.display()))?;
    }
    let existing = find_assembly_for_release(&args.datadir, &args.specie, &args.release)?;
    if existing.is_some() && !args.yes {
        bail!(
            "dataset already exists for {} release {}; pass --yes to rebuild",
            args.specie,
            args.release
        );
    }

    eprintln!(
        "Building dataset for {} release {} in {}",
        args.specie,
        args.release,
        args.datadir.display()
    );
    let dataset = ensembl::build_dataset(args)?;
    let path = save_dataset(args, &dataset)?;
    write_dataset_report(args, &dataset)?;
    if args.keep {
        write_transcriptome_fasta(args, &dataset)?;
    }
    println!("Dataset written: {}", path.display());
    Ok(())
}

pub fn update_dataset(args: &Args) -> Result<()> {
    let latest = ensembl::current_release(&args.specie)?;
    let mut updated = args.clone();
    updated.release = latest;
    if find_assembly_for_release(&updated.datadir, &updated.specie, &updated.release)?.is_some() {
        println!(
            "The last release for {} is {}, nothing to do.",
            updated.specie, updated.release
        );
        return Ok(());
    }
    make_dataset(&updated)
}

pub fn list_datasets(args: &Args) -> Result<()> {
    println!("Location of datasets: {}", args.datadir.display());
    if !args.datadir.exists() {
        println!("No datasets found");
        return Ok(());
    }
    let mut rows = Vec::new();
    for entry in std::fs::read_dir(&args.datadir)? {
        let entry = entry?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if let Some((species, assembly, release)) = parse_dataset_filename(&name) {
            rows.push((
                species.to_string(),
                assembly.to_string(),
                release.to_string(),
            ));
        }
    }
    rows.sort();
    if rows.is_empty() {
        println!("No datasets found");
    } else {
        println!("Datasets found:");
        for (species, assembly, release) in rows {
            println!("  - {species}: {release} ({assembly})");
        }
    }
    Ok(())
}

pub fn remove_dataset(args: &Args) -> Result<()> {
    let assembly = find_assembly_for_release(&args.datadir, &args.specie, &args.release)?
        .with_context(|| format!("dataset not found for {} {}", args.specie, args.release))?;
    let base = dataset_basename(&args.specie, &assembly, &args.release);
    for suffix in ["dataset.bin", "report.md", "transcriptome.fa"] {
        let path = args.datadir.join(format!("{base}.{suffix}"));
        if path.exists() {
            std::fs::remove_file(&path)
                .with_context(|| format!("failed to remove {}", path.display()))?;
            println!("removed {}", path.display());
        }
    }
    Ok(())
}

fn write_dataset_report(args: &Args, dataset: &Dataset) -> Result<()> {
    let base = dataset_basename(&dataset.species, &dataset.assembly, &dataset.release);
    let path = args.datadir.join(format!("{base}.report.md"));
    let text = format!(
        "# Kmerator-rs dataset\n\n- species: {}\n- release: {}\n- assembly: {}\n- genes: {}\n- transcripts: {}\n",
        dataset.species,
        dataset.release,
        dataset.assembly,
        dataset.geneinfo.genes.len(),
        dataset.transcriptome.records.len(),
    );
    std::fs::write(&path, text).with_context(|| format!("failed to write {}", path.display()))
}

fn write_transcriptome_fasta(args: &Args, dataset: &Dataset) -> Result<()> {
    let base = dataset_basename(&dataset.species, &dataset.assembly, &dataset.release);
    let path = args.datadir.join(format!("{base}.transcriptome.fa"));
    let file = std::fs::File::create(&path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = std::io::BufWriter::new(file);
    for rec in &dataset.transcriptome.records {
        crate::fastx::write_fasta_record(&mut writer, &rec.id, &rec.seq)?;
    }
    Ok(())
}

fn find_assembly_for_release(
    datadir: &Path,
    species: &str,
    release: &str,
) -> Result<Option<String>> {
    if !datadir.exists() {
        return Ok(None);
    }
    for entry in std::fs::read_dir(datadir)? {
        let entry = entry?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if let Some((sp, assembly, rel)) = parse_dataset_filename(&name)
            && sp == species
            && rel == release
        {
            return Ok(Some(assembly.to_string()));
        }
    }
    Ok(None)
}

fn parse_dataset_filename(name: &str) -> Option<(&str, &str, &str)> {
    let name = name.strip_suffix(".dataset.bin")?;
    let mut parts = name.split('.');
    let species = parts.next()?;
    let assembly = parts.next()?;
    let release = parts.next()?;
    if parts.next().is_some() {
        return None;
    }
    Some((species, assembly, release))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transcriptome_lookup_is_case_insensitive() {
        let transcriptome = Transcriptome::new(vec![TranscriptRecord {
            id: "ENST0001".to_string(),
            seq: b"ACGT".to_vec(),
        }]);

        assert_eq!(transcriptome.get("enst0001"), Some(&b"ACGT"[..]));
        assert_eq!(transcriptome.get("missing"), None);
    }

    #[test]
    fn dataset_basename_joins_species_assembly_and_release() {
        assert_eq!(
            dataset_basename("homo_sapiens", "GRCh38", "112"),
            "homo_sapiens.GRCh38.112"
        );
    }

    #[test]
    fn parse_dataset_filename_accepts_only_dataset_bins() {
        assert_eq!(
            parse_dataset_filename("homo_sapiens.GRCh38.112.dataset.bin"),
            Some(("homo_sapiens", "GRCh38", "112"))
        );
        assert_eq!(
            parse_dataset_filename("homo_sapiens.GRCh38.112.report.md"),
            None
        );
        assert_eq!(
            parse_dataset_filename("homo_sapiens.GRCh38.112.extra.dataset.bin"),
            None
        );
    }
}
