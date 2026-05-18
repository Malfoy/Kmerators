use crate::config::ConfigFile;
use anyhow::{Context, Result, bail};
use clap::Parser;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(
    name = "kmerators",
    version,
    about = "Find specific gene, transcript, or FASTA k-mers without Jellyfish",
    long_about = None
)]
pub struct Cli {
    /// Gene IDs, symbols, aliases, transcript IDs, or one file containing them.
    #[arg(short = 's', long = "selection", num_args = 1..)]
    pub selection: Vec<String>,

    /// FASTA/FASTQ file containing unannotated query sequences.
    #[arg(short = 'f', long = "fasta-file")]
    pub fasta_file: Option<PathBuf>,

    /// Storage directory for kmerator-rs datasets.
    #[arg(short = 'd', long = "datadir")]
    pub datadir: Option<PathBuf>,

    /// Reference genome FASTA/FASTQ to scan. Plain, gzip, zstd, and xz input is supported.
    #[arg(short = 'g', long = "genome")]
    pub genome: Option<PathBuf>,

    /// Local transcriptome FASTA/FASTQ for --fasta-file experiments, bypassing dataset lookup.
    #[arg(long = "transcriptome-fasta")]
    pub transcriptome_fasta: Option<PathBuf>,

    /// Ensembl species name or supported alias.
    #[arg(short = 'S', long = "specie")]
    pub specie: Option<String>,

    /// K-mer length. Can be repeated, for example: -k 31 -k 41 -k 51.
    #[arg(short = 'k', long = "kmer-length", action = clap::ArgAction::Append)]
    pub kmer_length: Vec<usize>,

    /// Minimizer length used for partitioning.
    #[arg(short = 'm', long = "minimizer-length")]
    pub minimizer_length: Option<usize>,

    /// Number of minimizer-routed hash tables.
    #[arg(long = "hash-tables", alias = "hash-table-count")]
    pub hash_table_count: Option<usize>,

    /// Ensembl release, or "last".
    #[arg(short = 'r', long = "release")]
    pub release: Option<String>,

    /// For genes, retain k-mers present in all isoforms.
    #[arg(long = "stringent")]
    pub stringent: bool,

    /// With --fasta-file, maximum count allowed in transcriptome.
    #[arg(short = 'T', long = "max-on-transcriptome")]
    pub max_on_transcriptome: Option<u32>,

    /// With --fasta-file, maximum count allowed in genome.
    #[arg(short = 'G', long = "max-on-genome")]
    pub max_on_genome: Option<u32>,

    /// Output directory.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Worker thread count.
    #[arg(short = 't', long = "thread")]
    pub thread: Option<usize>,

    /// Temporary directory.
    #[arg(long = "tmpdir")]
    pub tmpdir: Option<PathBuf>,

    /// Show more details.
    #[arg(short = 'D', long = "debug")]
    pub debug: bool,

    /// Keep intermediate files.
    #[arg(long = "keep")]
    pub keep: bool,

    /// Run non-interactively.
    #[arg(short = 'y', long = "yes")]
    pub yes: bool,

    /// Edit the configuration file.
    #[arg(short = 'e', long = "edit-config")]
    pub edit_config: bool,

    /// List local datasets.
    #[arg(short = 'l', long = "list-dataset", alias = "list-datasets")]
    pub list_dataset: bool,

    /// Remove a dataset.
    #[arg(long = "rm-dataset")]
    pub rm_dataset: bool,

    /// Build a dataset.
    #[arg(long = "mk-dataset")]
    pub mk_dataset: bool,

    /// Print the current Ensembl release.
    #[arg(long = "last-avail", alias = "last-available")]
    pub last_avail: bool,

    /// Build the latest dataset if absent.
    #[arg(short = 'u', long = "update-dataset")]
    pub update_dataset: bool,

    /// Print gene/transcript information.
    #[arg(long = "info", num_args = 1..)]
    pub info: Vec<String>,

    /// With --info, include transcript sequences.
    #[arg(short = 'a', long = "all")]
    pub all: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Extract,
    ListDataset,
    RemoveDataset,
    MakeDataset,
    LastAvailable,
    UpdateDataset,
    Info,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerSpec {
    pub kmer_length: usize,
    pub minimizer_length: usize,
}

#[derive(Debug, Clone)]
pub struct Args {
    pub mode: Mode,
    pub selection: Vec<String>,
    pub fasta_file: Option<PathBuf>,
    pub datadir: PathBuf,
    pub genome: Option<PathBuf>,
    pub transcriptome_fasta: Option<PathBuf>,
    pub specie: String,
    pub kmer_specs: Vec<KmerSpec>,
    pub hash_table_count: usize,
    pub release: String,
    pub stringent: bool,
    pub max_on_transcriptome: u32,
    pub max_on_genome: u32,
    pub output: PathBuf,
    pub thread: usize,
    pub tmpdir: Option<PathBuf>,
    pub debug: bool,
    pub keep: bool,
    pub yes: bool,
    pub info: Vec<String>,
    pub all: bool,
}

impl Args {
    pub fn from_cli_and_config(cli: Cli, cfg: &ConfigFile) -> Result<Self> {
        let mode = selected_mode(&cli)?;
        let kmer_lengths = kmer_lengths_from_cli_or_config(&cli, cfg)?;
        let explicit_minimizer_length = cli
            .minimizer_length
            .or_else(|| cfg.parse("minimizer_length"))
            .map(validate_minimizer_length)
            .transpose()?;
        let kmer_specs = build_kmer_specs(kmer_lengths, explicit_minimizer_length)?;
        let specie = cli
            .specie
            .or_else(|| cfg.get("specie"))
            .unwrap_or_else(|| "human".to_string());
        let specie = normalize_species(&specie);

        let release = cli
            .release
            .or_else(|| cfg.get("release"))
            .unwrap_or_else(|| "last".to_string());

        let datadir = cli
            .datadir
            .or_else(|| cfg.path("datadir"))
            .unwrap_or_else(|| PathBuf::from("./kmerators-data"));

        let hash_table_count = cli
            .hash_table_count
            .or_else(|| cfg.parse("hash_table_count"))
            .unwrap_or(1024);
        if hash_table_count == 0 {
            bail!("--hash-tables must be greater than zero");
        }

        let thread = cli
            .thread
            .or_else(|| cfg.parse("thread"))
            .unwrap_or_else(num_cpus::get)
            .max(1);

        let output = cli
            .output
            .or_else(|| cfg.path("output"))
            .unwrap_or_else(|| PathBuf::from("output"));

        let genome = cli.genome.or_else(|| cfg.path("genome"));
        let max_on_transcriptome = cli
            .max_on_transcriptome
            .or_else(|| cfg.parse("max_on_transcriptome"))
            .unwrap_or(0);
        let max_on_genome = cli
            .max_on_genome
            .or_else(|| cfg.parse("max_on_genome"))
            .unwrap_or(1);

        let stringent = cli.stringent || cfg.bool("stringent").unwrap_or(false);
        let keep = cli.keep || cfg.bool("keep").unwrap_or(false);
        let yes = cli.yes || cfg.bool("yes").unwrap_or(false);

        let mut args = Self {
            mode,
            selection: expand_selection(cli.selection)?,
            fasta_file: cli.fasta_file,
            datadir,
            genome,
            transcriptome_fasta: cli.transcriptome_fasta,
            specie,
            kmer_specs,
            hash_table_count,
            release,
            stringent,
            max_on_transcriptome,
            max_on_genome,
            output,
            thread,
            tmpdir: cli.tmpdir,
            debug: cli.debug,
            keep,
            yes,
            info: expand_selection(cli.info)?,
            all: cli.all,
        };

        validate(&mut args)?;
        Ok(args)
    }
}

fn kmer_lengths_from_cli_or_config(cli: &Cli, cfg: &ConfigFile) -> Result<Vec<usize>> {
    let mut lengths = if cli.kmer_length.is_empty() {
        vec![cfg.parse("kmer_length").unwrap_or(31)]
    } else {
        cli.kmer_length.clone()
    };

    let mut seen = HashSet::with_capacity(lengths.len());
    for &k in &lengths {
        if k == 0 {
            bail!("--kmer-length must be greater than zero");
        }
        if !seen.insert(k) {
            bail!("duplicate --kmer-length value: {k}");
        }
    }
    lengths.shrink_to_fit();
    Ok(lengths)
}

fn validate_minimizer_length(minimizer_length: usize) -> Result<usize> {
    if minimizer_length == 0 {
        bail!("--minimizer-length must be greater than zero");
    }
    if minimizer_length > 29 {
        bail!("--minimizer-length must be <= 29 for simd-minimizers");
    }
    Ok(minimizer_length)
}

fn build_kmer_specs(
    kmer_lengths: Vec<usize>,
    explicit_minimizer_length: Option<usize>,
) -> Result<Vec<KmerSpec>> {
    kmer_lengths
        .into_iter()
        .map(|kmer_length| {
            let minimizer_length = explicit_minimizer_length.unwrap_or_else(|| 9.min(kmer_length));
            if minimizer_length > kmer_length {
                bail!("--minimizer-length must be in 1..=k for every --kmer-length");
            }
            Ok(KmerSpec {
                kmer_length,
                minimizer_length,
            })
        })
        .collect()
}

fn selected_mode(cli: &Cli) -> Result<Mode> {
    let selected = [
        (
            !cli.selection.is_empty() || cli.fasta_file.is_some(),
            Mode::Extract,
        ),
        (cli.list_dataset, Mode::ListDataset),
        (cli.rm_dataset, Mode::RemoveDataset),
        (cli.mk_dataset, Mode::MakeDataset),
        (cli.last_avail, Mode::LastAvailable),
        (cli.update_dataset, Mode::UpdateDataset),
        (!cli.info.is_empty(), Mode::Info),
    ]
    .into_iter()
    .filter_map(|(on, mode)| on.then_some(mode))
    .collect::<Vec<_>>();

    match selected.as_slice() {
        [mode] => Ok(*mode),
        [] => {
            bail!("one action is required: --selection, --fasta-file, --info, or a dataset command")
        }
        _ => bail!("actions are mutually exclusive"),
    }
}

fn validate(args: &mut Args) -> Result<()> {
    if args.release == "last"
        && matches!(
            args.mode,
            Mode::Extract | Mode::Info | Mode::MakeDataset | Mode::UpdateDataset
        )
    {
        args.release = crate::ensembl::current_release(&args.specie).with_context(|| {
            format!(
                "failed to resolve current Ensembl release for {}",
                args.specie
            )
        })?;
    }

    match args.mode {
        Mode::Extract => {
            if args.selection.is_empty() == args.fasta_file.is_none() {
                bail!("use exactly one of --selection or --fasta-file");
            }
            let genome = args
                .genome
                .as_ref()
                .context("--genome is required for extraction")?;
            if !genome.is_file() {
                bail!("genome FASTA/FASTQ not found: {}", genome.display());
            }
            if let Some(transcriptome) = &args.transcriptome_fasta {
                if !args.selection.is_empty() {
                    bail!("--transcriptome-fasta currently supports --fasta-file extraction only");
                }
                if !transcriptome.is_file() {
                    bail!(
                        "transcriptome FASTA/FASTQ not found: {}",
                        transcriptome.display()
                    );
                }
            } else {
                ensure_datadir(&args.datadir, true)?;
            }
            if let Some(path) = &args.fasta_file
                && !path.is_file()
            {
                bail!("query FASTA/FASTQ not found: {}", path.display());
            }
        }
        Mode::Info => ensure_datadir(&args.datadir, true)?,
        Mode::ListDataset => ensure_datadir(&args.datadir, false)?,
        Mode::MakeDataset | Mode::UpdateDataset => ensure_datadir(&args.datadir, false)?,
        Mode::RemoveDataset => ensure_datadir(&args.datadir, true)?,
        Mode::LastAvailable => {}
    }

    Ok(())
}

fn ensure_datadir(path: &Path, must_exist: bool) -> Result<()> {
    if path.exists() {
        if !path.is_dir() {
            bail!("datadir is not a directory: {}", path.display());
        }
    } else if must_exist {
        bail!("datadir not found: {}", path.display());
    }
    Ok(())
}

fn expand_selection(values: Vec<String>) -> Result<Vec<String>> {
    if values.len() == 1 {
        let path = PathBuf::from(&values[0]);
        if path.is_file() {
            let text = std::fs::read_to_string(&path)
                .with_context(|| format!("failed to read {}", path.display()))?;
            let values = text
                .lines()
                .flat_map(|line| line.split('#').next().unwrap_or("").split_whitespace())
                .map(ToOwned::to_owned)
                .collect::<Vec<_>>();
            return Ok(values);
        }
    }
    Ok(values)
}

fn normalize_species(input: &str) -> String {
    match input.to_ascii_lowercase().as_str() {
        "human" => "homo_sapiens",
        "mouse" => "mus_musculus",
        "zebrafish" => "danio_rerio",
        "horse" => "equus_caballus",
        "hen" => "gallus_gallus",
        "c.elegans" => "caenorhabditis_elegans",
        "droso" => "drosophila_melanogaster",
        other => other,
    }
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;

    fn empty_cli() -> Cli {
        Cli {
            selection: Vec::new(),
            fasta_file: None,
            datadir: None,
            genome: None,
            transcriptome_fasta: None,
            specie: None,
            kmer_length: Vec::new(),
            minimizer_length: None,
            hash_table_count: None,
            release: None,
            stringent: false,
            max_on_transcriptome: None,
            max_on_genome: None,
            output: None,
            thread: None,
            tmpdir: None,
            debug: false,
            keep: false,
            yes: false,
            edit_config: false,
            list_dataset: false,
            rm_dataset: false,
            mk_dataset: false,
            last_avail: false,
            update_dataset: false,
            info: Vec::new(),
            all: false,
        }
    }

    #[test]
    fn fasta_with_local_transcriptome_does_not_require_datadir() {
        let tmp = tempfile::tempdir().unwrap();
        let query = tmp.path().join("query.fa");
        let transcriptome = tmp.path().join("transcriptome.fa");
        let genome = tmp.path().join("genome.fa");
        let output = tmp.path().join("out");
        std::fs::write(&query, b">q1\nACGTACCCC\n").unwrap();
        std::fs::write(&transcriptome, b">TR1\nGGGGTACCCAAAA\n").unwrap();
        std::fs::write(&genome, b">chr1\nTTTACGTAGGGACCCC\n").unwrap();

        let mut cli = empty_cli();
        cli.fasta_file = Some(query.clone());
        cli.transcriptome_fasta = Some(transcriptome.clone());
        cli.genome = Some(genome.clone());
        cli.specie = Some("toy_species".to_string());
        cli.release = Some("1".to_string());
        cli.kmer_length = vec![5];
        cli.minimizer_length = Some(3);
        cli.hash_table_count = Some(64);
        cli.output = Some(output.clone());
        cli.thread = Some(2);
        cli.yes = true;

        let cfg = config::path_for_tests(&tmp.path().join("config.ini"));
        let args = Args::from_cli_and_config(cli, &cfg).unwrap();

        assert_eq!(args.mode, Mode::Extract);
        assert_eq!(args.fasta_file, Some(query));
        assert_eq!(args.transcriptome_fasta, Some(transcriptome));
        assert_eq!(args.genome, Some(genome));
        assert_eq!(args.release, "1");
        assert_eq!(args.output, output);
        assert_eq!(args.kmer_specs.len(), 1);
        assert_eq!(args.kmer_specs[0].kmer_length, 5);
        assert_eq!(args.kmer_specs[0].minimizer_length, 3);
    }

    #[test]
    fn selection_file_expands_comments_and_whitespace() {
        let tmp = tempfile::tempdir().unwrap();
        let selection = tmp.path().join("selection.txt");
        std::fs::write(
            &selection,
            b"BRCA1  # symbol\n\nENST00000255409\tNPM1 # more queries\n",
        )
        .unwrap();

        let values = expand_selection(vec![selection.display().to_string()]).unwrap();

        assert_eq!(values, vec!["BRCA1", "ENST00000255409", "NPM1"]);
    }

    #[test]
    fn local_transcriptome_is_rejected_for_selection_mode() {
        let tmp = tempfile::tempdir().unwrap();
        let transcriptome = tmp.path().join("transcriptome.fa");
        let genome = tmp.path().join("genome.fa");
        std::fs::write(&transcriptome, b">TR1\nACGTAC\n").unwrap();
        std::fs::write(&genome, b">chr1\nACGTAC\n").unwrap();

        let mut cli = empty_cli();
        cli.selection = vec!["ENSG000001".to_string()];
        cli.transcriptome_fasta = Some(transcriptome);
        cli.genome = Some(genome);
        cli.release = Some("1".to_string());

        let cfg = config::path_for_tests(&tmp.path().join("config.ini"));
        let err = Args::from_cli_and_config(cli, &cfg).unwrap_err();

        assert!(
            err.to_string()
                .contains("--transcriptome-fasta currently supports --fasta-file extraction only")
        );
    }

    #[test]
    fn repeated_kmer_lengths_build_ordered_specs() {
        let mut cli = empty_cli();
        cli.list_dataset = true;
        cli.kmer_length = vec![31, 41, 51];

        let cfg = config::path_for_tests(Path::new("config.ini"));
        let args = Args::from_cli_and_config(cli, &cfg).unwrap();

        let lengths = args
            .kmer_specs
            .iter()
            .map(|spec| (spec.kmer_length, spec.minimizer_length))
            .collect::<Vec<_>>();
        assert_eq!(lengths, vec![(31, 9), (41, 9), (51, 9)]);
    }

    #[test]
    fn repeated_kmer_lengths_reject_duplicates() {
        let mut cli = empty_cli();
        cli.list_dataset = true;
        cli.kmer_length = vec![31, 31];

        let cfg = config::path_for_tests(Path::new("config.ini"));
        let err = Args::from_cli_and_config(cli, &cfg).unwrap_err();

        assert!(
            err.to_string()
                .contains("duplicate --kmer-length value: 31")
        );
    }

    #[test]
    fn species_aliases_are_normalized() {
        let mut cli = empty_cli();
        cli.list_dataset = true;
        cli.specie = Some("human".to_string());

        let cfg = config::path_for_tests(Path::new("config.ini"));
        let args = Args::from_cli_and_config(cli, &cfg).unwrap();

        assert_eq!(args.specie, "homo_sapiens");
    }

    #[test]
    fn mutually_exclusive_actions_are_rejected() {
        let mut cli = empty_cli();
        cli.selection = vec!["NPM1".to_string()];
        cli.list_dataset = true;

        let cfg = config::path_for_tests(Path::new("config.ini"));
        let err = Args::from_cli_and_config(cli, &cfg).unwrap_err();

        assert!(err.to_string().contains("actions are mutually exclusive"));
    }

    #[test]
    fn minimizer_length_must_not_exceed_kmer_length() {
        let mut cli = empty_cli();
        cli.list_dataset = true;
        cli.kmer_length = vec![3];
        cli.minimizer_length = Some(4);

        let cfg = config::path_for_tests(Path::new("config.ini"));
        let err = Args::from_cli_and_config(cli, &cfg).unwrap_err();

        assert!(
            err.to_string()
                .contains("--minimizer-length must be in 1..=k")
        );
    }
}
