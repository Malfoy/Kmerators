use crate::cli::{Args, KmerSpec, VERSION};
use crate::count::{self, CountIndex, CountKey, QuerySeed};
use crate::dataset::{self, Dataset, Gene, Transcriptome};
use crate::fastx;
use crate::kmer::{self, KeyMode};
use crate::output::{self, RunReport};
use anyhow::{Context, Result, bail};
use hashbrown::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ItemKind {
    Gene,
    Transcript,
    Fasta,
}

#[derive(Debug, Clone)]
struct QueryItem {
    given: String,
    id: String,
    ensg: Option<String>,
    enst: Option<String>,
    kind: ItemKind,
    seq: Vec<u8>,
    isoform_total: u32,
}

#[derive(Debug, Clone)]
struct Occurrence {
    item_idx: usize,
    pos: usize,
    seq: Vec<u8>,
    forward_key: u64,
    canonical_key: u64,
    transcriptome_count: u32,
    genome_count: u32,
    isoform_count: u32,
}

#[derive(Debug)]
struct SpecificRun {
    spec: KmerSpec,
    output: PathBuf,
    report: RunReport,
    occurrences: Vec<Occurrence>,
}

type FastaRecords = Vec<(String, Vec<u8>)>;

pub fn extract_specific_kmers(args: &Args) -> Result<()> {
    eprintln!("Load dataset {} / release {}", args.specie, args.release);
    if args.debug {
        let k_specs = args
            .kmer_specs
            .iter()
            .map(|spec| format!("k{}:m{}", spec.kmer_length, spec.minimizer_length))
            .collect::<Vec<_>>()
            .join(", ");
        eprintln!(
            "k/m={}, hash_tables={}, threads={}, tmpdir={}",
            k_specs,
            args.hash_table_count,
            args.thread,
            args.tmpdir
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "<auto>".to_string())
        );
    }
    let dataset = if let Some(path) = &args.transcriptome_fasta {
        load_local_transcriptome_dataset(args, path)?
    } else {
        dataset::load_dataset(args)?
    };
    let mut base_report = RunReport::default();

    let items = if !args.selection.is_empty() {
        resolve_selection(args, &dataset, &mut base_report)?
    } else {
        load_fasta_items(args, &mut base_report)?
    };
    if items.is_empty() {
        bail!("no usable query sequence");
    }

    eprintln!("Build query k-mer indexes");
    let mut runs = prepare_specific_runs(args, &items, &base_report)?;
    if runs.iter().all(|run| run.occurrences.is_empty()) {
        bail!("no query k-mers could be generated for requested k-mer lengths");
    }

    let seeds_by_run = runs
        .iter()
        .map(|run| seeds_for_occurrences(&run.occurrences))
        .collect::<Vec<_>>();

    eprintln!("Count query k-mers in transcriptome");
    let tr_indexes = runs
        .iter()
        .zip(&seeds_by_run)
        .map(|(run, seeds)| {
            Arc::new(CountIndex::new(
                run.spec.kmer_length,
                run.spec.minimizer_length,
                args.hash_table_count,
                CountKey::Forward,
                seeds,
            ))
        })
        .collect::<Vec<_>>();
    let tr_counts = count::count_transcriptome_records_many(
        &dataset.transcriptome.records,
        tr_indexes,
        args.thread,
    );

    eprintln!("Count query k-mers in genome");
    let genome = args.genome.as_ref().context("--genome is required")?;
    let ge_indexes = runs
        .iter()
        .zip(&seeds_by_run)
        .map(|(run, seeds)| {
            Arc::new(CountIndex::new(
                run.spec.kmer_length,
                run.spec.minimizer_length,
                args.hash_table_count,
                CountKey::Canonical,
                seeds,
            ))
        })
        .collect::<Vec<_>>();
    let ge_counts = count::count_path_many(genome, ge_indexes, args.thread)?;

    for ((run, tr_counts), ge_counts) in runs.iter_mut().zip(tr_counts).zip(ge_counts) {
        for (idx, occ) in run.occurrences.iter_mut().enumerate() {
            occ.transcriptome_count = tr_counts[idx];
            occ.genome_count = ge_counts[idx];
        }
    }

    for run in &mut runs {
        annotate_isoform_counts(run.spec, &dataset, &items, &mut run.occurrences)?;
    }

    eprintln!("Write output");
    std::fs::create_dir_all(&args.output)
        .with_context(|| format!("failed to create {}", args.output.display()))?;
    for run in &mut runs {
        write_run_outputs(args, &items, run)?;
    }

    let printed_report = if runs.len() == 1 {
        runs[0].report.clone()
    } else {
        let report = summarize_runs(&runs);
        output::write_markdown_report(
            &args.output.join("report.md"),
            &format!("kmerators v{VERSION}; k={}", format_k_list(&runs)),
            &args.specie,
            &args.release,
            &report,
        )?;
        report
    };
    output::print_report(&printed_report);
    Ok(())
}

fn load_local_transcriptome_dataset(args: &Args, path: &std::path::Path) -> Result<Dataset> {
    let records = fastx::read_records(path)?
        .into_iter()
        .map(|rec| dataset::TranscriptRecord {
            id: fastx::first_token_without_version(&rec.header),
            seq: rec.seq,
        })
        .collect::<Vec<_>>();
    Ok(Dataset {
        species: args.specie.clone(),
        release: args.release.clone(),
        assembly: "local".to_string(),
        geneinfo: dataset::GeneInfo {
            assembly: "local".to_string(),
            chromosomes: Vec::new(),
            version: dataset::DATASET_VERSION,
            genes: std::collections::HashMap::new(),
            symbols: std::collections::HashMap::new(),
            aliases: std::collections::HashMap::new(),
            transcripts: std::collections::HashMap::new(),
        },
        transcriptome: dataset::Transcriptome::new(records),
    })
}

fn prepare_specific_runs(
    args: &Args,
    items: &[QueryItem],
    base_report: &RunReport,
) -> Result<Vec<SpecificRun>> {
    let multi_k = args.kmer_specs.len() > 1;
    args.kmer_specs
        .iter()
        .map(|&spec| {
            let mut report = base_report.clone();
            let occurrences = build_occurrences(spec, items, &mut report)?;
            Ok(SpecificRun {
                spec,
                output: output_dir_for_spec(args, spec, multi_k),
                report,
                occurrences,
            })
        })
        .collect()
}

fn output_dir_for_spec(args: &Args, spec: KmerSpec, multi_k: bool) -> PathBuf {
    if multi_k {
        args.output.join(format!("k{}", spec.kmer_length))
    } else {
        args.output.clone()
    }
}

fn seeds_for_occurrences(occurrences: &[Occurrence]) -> Vec<QuerySeed> {
    occurrences
        .iter()
        .enumerate()
        .map(|(id, occ)| QuerySeed {
            id,
            seq: occ.seq.clone(),
            forward_key: occ.forward_key,
            canonical_key: occ.canonical_key,
        })
        .collect()
}

fn write_run_outputs(args: &Args, items: &[QueryItem], run: &mut SpecificRun) -> Result<()> {
    let (kmers, contigs, masked) = build_outputs(args, items, &run.occurrences, &mut run.report);
    std::fs::create_dir_all(&run.output)
        .with_context(|| format!("failed to create {}", run.output.display()))?;
    output::write_fasta(&run.output.join("kmers.fa"), &kmers)?;
    output::write_fasta(&run.output.join("contigs.fa"), &contigs)?;
    output::write_fasta(&run.output.join("masked.fa"), &masked)?;
    let command = if args.kmer_specs.len() == 1 {
        format!("kmerators v{VERSION}")
    } else {
        format!("kmerators v{VERSION}; k={}", run.spec.kmer_length)
    };
    output::write_markdown_report(
        &run.output.join("report.md"),
        &command,
        &args.specie,
        &args.release,
        &run.report,
    )
}

fn summarize_runs(runs: &[SpecificRun]) -> RunReport {
    let mut report = RunReport::default();
    for run in runs {
        let prefix = format!("k={}: ", run.spec.kmer_length);
        append_prefixed(&mut report.done, &run.report.done, &prefix);
        append_prefixed(&mut report.multiple, &run.report.multiple, &prefix);
        append_prefixed(&mut report.failed, &run.report.failed, &prefix);
        append_prefixed(&mut report.warning, &run.report.warning, &prefix);
    }
    report
}

fn append_prefixed(dst: &mut Vec<String>, src: &[String], prefix: &str) {
    dst.extend(src.iter().map(|line| format!("{prefix}{line}")));
}

fn format_k_list(runs: &[SpecificRun]) -> String {
    runs.iter()
        .map(|run| run.spec.kmer_length.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

pub fn show_info(args: &Args) -> Result<()> {
    let dataset = dataset::load_dataset(args)?;
    for query in &args.info {
        let upper = query.to_ascii_uppercase();
        if let Some(ensgs) = dataset.geneinfo.symbols.get(&upper) {
            print_gene_group(args, &dataset, query, "gene symbol", ensgs)?;
        } else if let Some(ensgs) = dataset.geneinfo.aliases.get(&upper) {
            print_gene_group(args, &dataset, query, "alias", ensgs)?;
        } else if let Some(gene) = dataset.geneinfo.genes.get(&upper) {
            print_gene(args, &dataset, query, "ensembl gene name", &upper, gene)?;
        } else if let Some(ensg) = dataset.geneinfo.transcripts.get(&upper) {
            let gene = dataset
                .geneinfo
                .genes
                .get(ensg)
                .with_context(|| format!("parent gene {ensg} not found"))?;
            print_gene(args, &dataset, query, "transcript", ensg, gene)?;
        } else {
            println!("\n=== {query} ===\n  not found");
        }
    }
    Ok(())
}

fn resolve_selection(
    args: &Args,
    dataset: &Dataset,
    report: &mut RunReport,
) -> Result<Vec<QueryItem>> {
    let mut items = Vec::new();
    for given in &args.selection {
        if given.starts_with("ENS") && given.contains('.') {
            report.failed.push(format!(
                "{given}: Ensembl IDs with version suffix are not supported"
            ));
            continue;
        }
        let upper = given.to_ascii_uppercase();
        let mut matches = Vec::<String>::new();
        let mut source = ItemKind::Gene;

        if let Some(ensgs) = dataset.geneinfo.symbols.get(&upper) {
            matches.extend(ensgs.clone());
        } else if let Some(ensgs) = dataset.geneinfo.aliases.get(&upper) {
            matches.extend(ensgs.clone());
        } else if dataset.geneinfo.genes.contains_key(&upper) {
            matches.push(upper.clone());
        } else if let Some(ensg) = dataset.geneinfo.transcripts.get(&upper) {
            source = ItemKind::Transcript;
            matches.push(ensg.clone());
        }

        if matches.is_empty() {
            report
                .failed
                .push(format!("{given}: not found in transcriptome"));
            continue;
        }
        if matches.len() > 1 {
            report.multiple.push(format!(
                "{given}: {} ({})",
                matches.len(),
                matches.join(", ")
            ));
        }

        for ensg in matches {
            let gene = dataset
                .geneinfo
                .genes
                .get(&ensg)
                .with_context(|| format!("gene metadata missing for {ensg}"))?;
            let enst = if source == ItemKind::Transcript {
                upper.clone()
            } else {
                match &gene.canonical {
                    Some(enst) => enst.clone(),
                    None => {
                        report
                            .failed
                            .push(format!("{given}: no canonical transcript for {ensg}"));
                        continue;
                    }
                }
            };
            let Some(seq) = dataset.transcriptome.get(&enst) else {
                report.failed.push(format!(
                    "{given}: transcript not found in transcriptome ({enst})"
                ));
                continue;
            };
            items.push(QueryItem {
                given: given.clone(),
                id: if source == ItemKind::Transcript {
                    enst.clone()
                } else {
                    format!("{}:{}", gene.symbol.as_deref().unwrap_or(&ensg), enst)
                },
                ensg: Some(ensg.clone()),
                enst: Some(enst),
                kind: source,
                seq: seq.to_vec(),
                isoform_total: gene.transcripts.len() as u32,
            });
        }
    }
    Ok(items)
}

fn load_fasta_items(args: &Args, _report: &mut RunReport) -> Result<Vec<QueryItem>> {
    let path = args
        .fasta_file
        .as_ref()
        .context("--fasta-file is required")?;
    let mut items = Vec::new();
    for rec in fastx::read_records(path)? {
        let id = fastx::first_token(&rec.header);
        items.push(QueryItem {
            given: id.clone(),
            id,
            ensg: None,
            enst: None,
            kind: ItemKind::Fasta,
            seq: rec.seq,
            isoform_total: 0,
        });
    }
    Ok(items)
}

fn build_occurrences(
    spec: KmerSpec,
    items: &[QueryItem],
    report: &mut RunReport,
) -> Result<Vec<Occurrence>> {
    let kmer_length = spec.kmer_length;
    let mode = KeyMode::for_k(kmer_length);
    let mut occurrences = Vec::new();
    for (item_idx, item) in items.iter().enumerate() {
        let before = occurrences.len();
        for (chunk_start, chunk) in kmer::valid_chunks(&item.seq) {
            if chunk.len() < kmer_length {
                continue;
            }
            let keys = kmer::kmer_keys_for_chunk(chunk, kmer_length, mode);
            for (offset, (forward_key, canonical_key)) in keys.into_iter().enumerate() {
                let pos = chunk_start + offset + 1;
                occurrences.push(Occurrence {
                    item_idx,
                    pos,
                    seq: chunk[offset..offset + kmer_length].to_vec(),
                    forward_key,
                    canonical_key,
                    transcriptome_count: 0,
                    genome_count: 0,
                    isoform_count: 0,
                });
            }
        }
        if before == occurrences.len() {
            if item.seq.len() < kmer_length {
                report.failed.push(format!(
                    "{}: sequence too short ({} < {})",
                    item.given,
                    item.seq.len(),
                    kmer_length
                ));
            } else {
                report.failed.push(format!(
                    "{}: no valid ACTG k-mer generated for k={}",
                    item.given, kmer_length
                ));
            }
        }
    }
    Ok(occurrences)
}

fn annotate_isoform_counts(
    spec: KmerSpec,
    dataset: &Dataset,
    items: &[QueryItem],
    occurrences: &mut [Occurrence],
) -> Result<()> {
    let mut cache = HashMap::<String, HashMap<u64, u32>>::new();
    for occ in occurrences {
        let item = &items[occ.item_idx];
        if item.kind != ItemKind::Gene {
            continue;
        }
        let ensg = item.ensg.as_ref().expect("gene item has ensg");
        if !cache.contains_key(ensg) {
            let gene = dataset
                .geneinfo
                .genes
                .get(ensg)
                .with_context(|| format!("gene metadata missing for {ensg}"))?;
            let counts = isoform_presence_counts(spec, &dataset.transcriptome, gene);
            cache.insert(ensg.clone(), counts);
        }
        occ.isoform_count = cache
            .get(ensg)
            .and_then(|counts| counts.get(&occ.forward_key))
            .copied()
            .unwrap_or(0);
    }
    Ok(())
}

fn isoform_presence_counts(
    spec: KmerSpec,
    transcriptome: &Transcriptome,
    gene: &Gene,
) -> HashMap<u64, u32> {
    let kmer_length = spec.kmer_length;
    let mode = KeyMode::for_k(kmer_length);
    let mut counts = HashMap::<u64, u32>::new();
    for enst in &gene.transcripts {
        let Some(seq) = transcriptome.get(enst) else {
            continue;
        };
        let mut seen = HashSet::<u64>::new();
        for (_, chunk) in kmer::valid_chunks(seq) {
            if chunk.len() < kmer_length {
                continue;
            }
            for (forward, _) in kmer::kmer_keys_for_chunk(chunk, kmer_length, mode) {
                seen.insert(forward);
            }
        }
        for key in seen {
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    counts
}

fn build_outputs(
    args: &Args,
    items: &[QueryItem],
    occurrences: &[Occurrence],
    report: &mut RunReport,
) -> (FastaRecords, FastaRecords, FastaRecords) {
    let mut kmers = Vec::new();
    let mut contigs = Vec::new();
    let mut masked = Vec::new();

    for (item_idx, item) in items.iter().enumerate() {
        let item_occs = occurrences
            .iter()
            .filter(|occ| occ.item_idx == item_idx)
            .collect::<Vec<_>>();
        if item_occs.is_empty() {
            continue;
        }
        let mut specific_count = 0usize;
        let mut contig_count = 0usize;
        let mut current_contig = Vec::<u8>::new();
        let mut contig_start = 0usize;
        let mut current_contig_no = 0usize;
        let mut last_pos = 0usize;

        for occ in item_occs {
            if is_specific(args, item, occ) {
                specific_count += 1;
                if current_contig.is_empty() || occ.pos != last_pos + 1 {
                    if !current_contig.is_empty() {
                        contigs.push((
                            contig_header(item, current_contig_no, contig_start),
                            std::mem::take(&mut current_contig),
                        ));
                    }
                    contig_count += 1;
                    current_contig_no = contig_count;
                    contig_start = occ.pos;
                    current_contig = occ.seq.clone();
                } else {
                    current_contig.push(*occ.seq.last().expect("k-mer is non-empty"));
                }
                last_pos = occ.pos;
                kmers.push((kmer_header(item, occ, current_contig_no), occ.seq.clone()));
            } else {
                masked.push((masked_header(item, occ), occ.seq.clone()));
            }
        }

        if !current_contig.is_empty() {
            contigs.push((
                contig_header(item, current_contig_no, contig_start),
                current_contig,
            ));
        }

        if specific_count == 0 {
            report
                .failed
                .push(format!("{}: no specific kmers found", item.given));
        } else {
            report.done.push(format!(
                "{}: {} - kmers/contigs: {}/{} ({})",
                item.given,
                item.id,
                specific_count,
                contig_count,
                item_kind_label(item.kind)
            ));
        }
    }

    (kmers, contigs, masked)
}

fn is_specific(args: &Args, item: &QueryItem, occ: &Occurrence) -> bool {
    match item.kind {
        ItemKind::Gene => {
            if occ.genome_count > 1 {
                return false;
            }
            if args.stringent {
                occ.transcriptome_count == item.isoform_total
                    && occ.isoform_count == item.isoform_total
            } else {
                occ.transcriptome_count > 0 && occ.transcriptome_count == occ.isoform_count
            }
        }
        ItemKind::Transcript => occ.transcriptome_count == 1 && occ.genome_count <= 1,
        ItemKind::Fasta => {
            occ.transcriptome_count <= args.max_on_transcriptome
                && occ.genome_count <= args.max_on_genome
        }
    }
}

fn kmer_header(item: &QueryItem, occ: &Occurrence, contig_no: usize) -> String {
    match item.kind {
        ItemKind::Gene => format!(
            "{}:{}.kmer{} ct:{} tr:{}/{}",
            item.given.to_ascii_uppercase(),
            item.enst.as_deref().unwrap_or("NA"),
            occ.pos,
            contig_no,
            occ.isoform_count,
            item.isoform_total
        ),
        ItemKind::Transcript => format!(
            "{}:{}.kmer{} ct:{}",
            item.given.to_ascii_uppercase(),
            item.enst.as_deref().unwrap_or(&item.id),
            occ.pos,
            contig_no
        ),
        ItemKind::Fasta => format!("{}.kmer{} ct:{}", item.id, occ.pos, contig_no),
    }
}

fn contig_header(item: &QueryItem, contig_no: usize, pos: usize) -> String {
    match item.kind {
        ItemKind::Gene => format!(
            "{}:{}.contig_{} (at position {})",
            item.given.to_ascii_uppercase(),
            item.enst.as_deref().unwrap_or("NA"),
            contig_no,
            pos
        ),
        ItemKind::Transcript => format!(
            "{}.contig_{} (at position {})",
            item.enst.as_deref().unwrap_or(&item.id),
            contig_no,
            pos
        ),
        ItemKind::Fasta => format!("{}.contig_{} (at position {})", item.id, contig_no, pos),
    }
}

fn masked_header(item: &QueryItem, occ: &Occurrence) -> String {
    match item.kind {
        ItemKind::Gene => format!(
            "{}:{}.kmer{} tr:{}/{} genome:{} transcriptome:{}",
            item.given.to_ascii_uppercase(),
            item.enst.as_deref().unwrap_or("NA"),
            occ.pos,
            occ.isoform_count,
            item.isoform_total,
            occ.genome_count,
            occ.transcriptome_count
        ),
        ItemKind::Transcript => format!(
            "{}:{}.kmer{} genome:{} transcriptome:{}",
            item.given.to_ascii_uppercase(),
            item.enst.as_deref().unwrap_or(&item.id),
            occ.pos,
            occ.genome_count,
            occ.transcriptome_count
        ),
        ItemKind::Fasta => format!(
            "{}.kmer{} genome:{} transcriptome:{}",
            item.id, occ.pos, occ.genome_count, occ.transcriptome_count
        ),
    }
}

fn item_kind_label(kind: ItemKind) -> &'static str {
    match kind {
        ItemKind::Gene => "gene",
        ItemKind::Transcript => "transcript",
        ItemKind::Fasta => "fasta",
    }
}

fn print_gene_group(
    args: &Args,
    dataset: &Dataset,
    query: &str,
    kind: &str,
    ensgs: &[String],
) -> Result<()> {
    println!("\n=== {query} ({kind}) ({} found) ===", ensgs.len());
    for ensg in ensgs {
        let gene = dataset
            .geneinfo
            .genes
            .get(ensg)
            .with_context(|| format!("gene metadata missing for {ensg}"))?;
        print_gene(args, dataset, query, kind, ensg, gene)?;
    }
    Ok(())
}

fn print_gene(
    args: &Args,
    dataset: &Dataset,
    query: &str,
    kind: &str,
    ensg: &str,
    gene: &Gene,
) -> Result<()> {
    println!("\n=== {query} ({kind}) ===");
    println!("  Ensembl ID           {ensg}");
    println!(
        "  Gene Name            {}",
        gene.symbol.as_deref().unwrap_or("")
    );
    println!("  Specie               {}", args.specie);
    println!("  Assembly             {}", dataset.assembly);
    println!(
        "  Coordinates          {}:{}-{}",
        gene.chr, gene.start, gene.end
    );
    println!("  Strand               {}", gene.strand);
    println!(
        "  Canonical transcript {}",
        gene.canonical.as_deref().unwrap_or("unknown")
    );
    println!("  Biotype              {}", gene.biotype);
    println!("  Description          {}", gene.desc);
    println!(
        "  Transcripts ({})     {}",
        gene.transcripts.len(),
        gene.transcripts.join(" ")
    );
    if args.all {
        for enst in &gene.transcripts {
            if let Some(seq) = dataset.transcriptome.get(enst) {
                println!("{enst} ({})", seq.len());
                println!("{}", String::from_utf8_lossy(seq));
            }
        }
    }
    Ok(())
}

#[allow(dead_code)]
fn output_path(args: &Args, name: &str) -> PathBuf {
    args.output.join(name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Mode;
    use crate::dataset::{GeneInfo, TranscriptRecord};
    use std::collections::HashMap as StdHashMap;

    #[test]
    fn extracts_specific_kmers_from_fasta_query() {
        let tmp = tempfile::tempdir().unwrap();
        let datadir = tmp.path().join("data");
        let output = tmp.path().join("out");
        let genome = tmp.path().join("genome.fa");
        let query = tmp.path().join("query.fa");
        std::fs::create_dir_all(&datadir).unwrap();
        std::fs::write(&genome, b">chr1\nACGTACGT\n").unwrap();
        std::fs::write(&query, b">q1\nACGTAC\n").unwrap();

        let args = Args {
            mode: Mode::Extract,
            selection: Vec::new(),
            fasta_file: Some(query),
            datadir: datadir.clone(),
            genome: Some(genome),
            transcriptome_fasta: None,
            specie: "test_species".to_string(),
            kmer_specs: vec![KmerSpec {
                kmer_length: 3,
                minimizer_length: 2,
            }],
            hash_table_count: 1024,
            release: "1".to_string(),
            stringent: false,
            max_on_transcriptome: 0,
            max_on_genome: 10,
            output: output.clone(),
            thread: 2,
            tmpdir: None,
            debug: false,
            keep: false,
            yes: true,
            info: Vec::new(),
            all: false,
        };

        let dataset = Dataset {
            species: args.specie.clone(),
            release: args.release.clone(),
            assembly: "ASM".to_string(),
            geneinfo: GeneInfo {
                assembly: "ASM".to_string(),
                chromosomes: vec!["chr1".to_string()],
                version: dataset::DATASET_VERSION,
                genes: StdHashMap::new(),
                symbols: StdHashMap::new(),
                aliases: StdHashMap::new(),
                transcripts: StdHashMap::new(),
            },
            transcriptome: dataset::Transcriptome::new(vec![TranscriptRecord {
                id: "TR1".to_string(),
                seq: b"TTTTTT".to_vec(),
            }]),
        };
        dataset::save_dataset(&args, &dataset).unwrap();

        extract_specific_kmers(&args).unwrap();

        let kmers = std::fs::read_to_string(output.join("kmers.fa")).unwrap();
        let contigs = std::fs::read_to_string(output.join("contigs.fa")).unwrap();
        assert!(kmers.contains(">q1.kmer1 ct:1"));
        assert!(contigs.contains(">q1.contig_1 (at position 1)\nACGTAC"));
    }
}
