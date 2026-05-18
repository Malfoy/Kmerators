use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::FromStr;

const DEFAULT_CONFIG: &str = r#"[CMD_ARGS]
## --datadir option
# datadir = /path/to/kmerators/directory

## --release option: release number or last
# release = last

## --genome option: path to reference genome FASTA/FASTQ
# genome = /genomes/GRCh38.fa.gz

## --specie option
# specie = human

## --thread option
# thread = 1

## --kmer-length option
# kmer_length = 31

## --minimizer-length option
# minimizer_length = 9

## number of minimizer-routed hash tables
# hash_table_count = 1024

## --max-on-transcriptome, only with --fasta-file
# max_on_transcriptome = 0

## --max-on-genome, only with --fasta-file
# max_on_genome = 1

## --stringent option, only with --selection
# stringent = false

## output directory
# output = ./output

## assumes yes as prompt answer
# yes = false

## keep intermediate files
# keep = false
"#;

#[derive(Debug, Clone)]
pub struct ConfigFile {
    path: PathBuf,
    values: HashMap<String, String>,
}

impl ConfigFile {
    pub fn load_or_create(appname: &str) -> Result<Self> {
        let path = config_path(appname)?;
        if !path.exists() {
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }
            let mut file = std::fs::File::create(&path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            file.write_all(DEFAULT_CONFIG.as_bytes())
                .with_context(|| format!("failed to write {}", path.display()))?;
        }

        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        Ok(Self {
            path,
            values: parse_cmd_args_section(&text),
        })
    }

    pub fn get(&self, key: &str) -> Option<String> {
        self.values.get(&normalize_key(key)).cloned()
    }

    pub fn path(&self, key: &str) -> Option<PathBuf> {
        self.get(key).map(PathBuf::from)
    }

    pub fn parse<T>(&self, key: &str) -> Option<T>
    where
        T: FromStr,
    {
        self.get(key).and_then(|value| value.parse::<T>().ok())
    }

    pub fn bool(&self, key: &str) -> Option<bool> {
        self.get(key)
            .and_then(|value| match value.to_ascii_lowercase().as_str() {
                "true" | "yes" | "1" | "on" => Some(true),
                "false" | "no" | "0" | "off" => Some(false),
                _ => None,
            })
    }

    pub fn edit(&self) -> Result<()> {
        let editor = std::env::var("EDITOR")
            .ok()
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "editor".to_string());
        Command::new(editor)
            .arg(&self.path)
            .status()
            .with_context(|| format!("failed to run editor for {}", self.path.display()))?;
        Ok(())
    }
}

fn config_path(appname: &str) -> Result<PathBuf> {
    let dir = if unsafe { libc_geteuid() } == 0 {
        PathBuf::from("/etc").join(appname)
    } else {
        dirs::config_dir()
            .context("could not determine user configuration directory")?
            .join(appname)
    };
    Ok(dir.join("config-v3.ini"))
}

unsafe fn libc_geteuid() -> u32 {
    unsafe extern "C" {
        fn geteuid() -> u32;
    }
    unsafe { geteuid() }
}

fn parse_cmd_args_section(text: &str) -> HashMap<String, String> {
    let mut in_section = false;
    let mut values = HashMap::new();

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with(';') {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            in_section = &line[1..line.len() - 1] == "CMD_ARGS";
            continue;
        }
        if !in_section {
            continue;
        }
        if let Some((key, value)) = line.split_once('=') {
            values.insert(normalize_key(key.trim()), value.trim().to_string());
        }
    }

    values
}

fn normalize_key(key: &str) -> String {
    key.trim_start_matches('-').replace('-', "_")
}

#[allow(dead_code)]
pub fn path_for_tests(path: &Path) -> ConfigFile {
    ConfigFile {
        path: path.to_path_buf(),
        values: HashMap::new(),
    }
}
