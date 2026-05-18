mod cli;
mod config;
mod count;
mod dataset;
mod ensembl;
mod fastx;
mod kmer;
mod output;
mod specific;

use anyhow::{Context, Result};
use clap::Parser;

fn main() -> Result<()> {
    let raw = cli::Cli::parse();
    let cfg = config::ConfigFile::load_or_create("kmerators")
        .context("failed to load kmerators configuration")?;

    if raw.edit_config {
        cfg.edit().context("failed to open configuration file")?;
        return Ok(());
    }

    let args = cli::Args::from_cli_and_config(raw, &cfg)?;

    match args.mode {
        cli::Mode::ListDataset => dataset::list_datasets(&args),
        cli::Mode::RemoveDataset => dataset::remove_dataset(&args),
        cli::Mode::MakeDataset => dataset::make_dataset(&args),
        cli::Mode::LastAvailable => {
            let release = ensembl::current_release(&args.specie)?;
            println!("{release}");
            Ok(())
        }
        cli::Mode::UpdateDataset => dataset::update_dataset(&args),
        cli::Mode::Info => specific::show_info(&args),
        cli::Mode::Extract => specific::extract_specific_kmers(&args),
    }
}
