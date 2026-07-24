use crate::cli::{Args, KmerSpec};
use crate::count::{self, CountIndex, CountKey, QuerySeed};
use crate::dataset::{self, Dataset, Gene, Transcriptome};
use crate::fastx;
use crate::kmer::{self, KeyMode};
use crate::output::{self, RunReport};
use anyhow::{Context, Result, bail};
use hashbrown::{HashMap, HashSet};
use std::ffi::OsStr;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, Instant, UNIX_EPOCH};

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

pub fn extract_specific_kmers(args: &Args) -> Result<()> {
    let total_timer = PhaseTimer::start("Total extraction");
    let mut phase_timings = Vec::new();
    if args.debug {
        let k_specs = args
            .kmer_specs
            .iter()
            .map(|spec| format!("k{}:m{}", spec.kmer_length, spec.minimizer_length))
            .collect::<Vec<_>>()
            .join(", ");
        eprintln!(
            "k/m={}, hash_tables={}, threads={}",
            k_specs, args.hash_table_count, args.thread
        );
    }

    let dataset = {
        let timer = PhaseTimer::start(format!(
            "Load dataset {} / release {}",
            args.specie, args.release
        ));
        let dataset = if args.transcriptome_fasta.is_some() {
            local_dataset(args)
        } else {
            dataset::load_dataset(args)?
        };
        phase_timings.push(timer.finish());
        dataset
    };
    let mut base_report = RunReport::default();

    let items = {
        let timer = PhaseTimer::start("Resolve query sequences");
        let items = if !args.selection.is_empty() {
            resolve_selection(args, &dataset, &mut base_report)?
        } else {
            load_fasta_items(args, &mut base_report)?
        };
        phase_timings.push(timer.finish());
        items
    };
    if items.is_empty() {
        bail!("no usable query sequence");
    }

    let mut runs = {
        let timer = PhaseTimer::start("Build query k-mer occurrences");
        let runs = prepare_specific_runs(args, &items, &base_report)?;
        phase_timings.push(timer.finish());
        runs
    };
    if runs.iter().all(|run| run.occurrences.is_empty()) {
        bail!("no query k-mers could be generated for requested k-mer lengths");
    }

    let tr_indexes = {
        let timer = PhaseTimer::start("Build transcriptome count indexes");
        let indexes = runs
            .iter()
            .map(|run| {
                let seeds = seeds_for_occurrences(run.spec, &items, &run.occurrences);
                Arc::new(CountIndex::new(
                    run.spec.kmer_length,
                    run.spec.minimizer_length,
                    args.hash_table_count,
                    CountKey::Forward,
                    &seeds,
                ))
            })
            .collect::<Vec<_>>();
        phase_timings.push(timer.finish());
        indexes
    };
    let tr_counts = {
        let timer = PhaseTimer::start("Count query k-mers in transcriptome");
        let counts = if let Some(path) = &args.transcriptome_fasta {
            count::count_path_many(path, tr_indexes, args.thread)?
        } else {
            count::count_transcriptome_records_many(
                &dataset.transcriptome.records,
                tr_indexes,
                args.thread,
            )
        };
        phase_timings.push(timer.finish());
        counts
    };

    let genome = args.genome.as_ref().context("--genome is required")?;
    let ge_indexes = {
        let timer = PhaseTimer::start("Build genome count indexes");
        let indexes = runs
            .iter()
            .map(|run| {
                let seeds = seeds_for_occurrences(run.spec, &items, &run.occurrences);
                Arc::new(CountIndex::new(
                    run.spec.kmer_length,
                    run.spec.minimizer_length,
                    args.hash_table_count,
                    CountKey::Canonical,
                    &seeds,
                ))
            })
            .collect::<Vec<_>>();
        phase_timings.push(timer.finish());
        indexes
    };
    let ge_counts = {
        let timer = PhaseTimer::start("Count query k-mers in genome");
        let counts = count::count_path_many(genome, ge_indexes, args.thread)?;
        phase_timings.push(timer.finish());
        counts
    };

    {
        let timer = PhaseTimer::start("Assign k-mer counts");
        for ((run, tr_counts), ge_counts) in runs.iter_mut().zip(tr_counts).zip(ge_counts) {
            for (idx, occ) in run.occurrences.iter_mut().enumerate() {
                occ.transcriptome_count = tr_counts[idx];
                occ.genome_count = ge_counts[idx];
            }
        }
        phase_timings.push(timer.finish());
    }

    {
        let timer = PhaseTimer::start("Annotate isoform counts");
        for run in &mut runs {
            annotate_isoform_counts(run.spec, &dataset, &items, &mut run.occurrences)?;
        }
        phase_timings.push(timer.finish());
    }

    {
        let timer = PhaseTimer::start("Write output files");
        std::fs::create_dir_all(&args.output)
            .with_context(|| format!("failed to create {}", args.output.display()))?;
        for run in &mut runs {
            let output_timer =
                PhaseTimer::start(format!("Write output for k={}", run.spec.kmer_length));
            write_run_outputs(args, &items, run)?;
            phase_timings.push(output_timer.finish());
        }
        phase_timings.push(timer.finish());
    }

    let printed_report = {
        let timer = PhaseTimer::start("Write summary report");
        let inputs = collect_input_metadata(args)?;
        let total_elapsed = total_timer.elapsed();
        for run in &runs {
            let metadata = report_metadata(
                args,
                &inputs,
                &phase_timings,
                &[(run.spec.kmer_length, run.spec.minimizer_length)],
                total_elapsed,
            );
            output::write_markdown_report(&run.output.join("report.md"), &metadata, &run.report)?;
        }
        let printed_report = if runs.len() == 1 {
            runs[0].report.clone()
        } else {
            let report = summarize_runs(&runs);
            let specs = runs
                .iter()
                .map(|run| (run.spec.kmer_length, run.spec.minimizer_length))
                .collect::<Vec<_>>();
            let metadata = report_metadata(args, &inputs, &phase_timings, &specs, total_elapsed);
            output::write_markdown_report(&args.output.join("report.md"), &metadata, &report)?;
            report
        };
        phase_timings.push(timer.finish());
        printed_report
    };
    output::print_report(&printed_report);
    total_timer.finish();
    Ok(())
}

#[derive(Debug)]
struct PhaseTimer {
    label: String,
    start: Instant,
    finished: bool,
}

impl PhaseTimer {
    fn start(label: impl Into<String>) -> Self {
        let label = label.into();
        eprintln!("{label}...");
        Self {
            label,
            start: Instant::now(),
            finished: false,
        }
    }

    fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }

    fn finish(mut self) -> output::PhaseTiming {
        self.finished = true;
        let elapsed = self.start.elapsed();
        eprintln!("{} done in {}", self.label, format_duration(elapsed));
        output::PhaseTiming {
            label: self.label.clone(),
            elapsed,
        }
    }
}

impl Drop for PhaseTimer {
    fn drop(&mut self) {
        if !self.finished {
            eprintln!(
                "{} failed after {}",
                self.label,
                format_duration(self.start.elapsed())
            );
        }
    }
}

fn format_duration(duration: Duration) -> String {
    format!("{:.3}s", duration.as_secs_f64())
}

fn collect_input_metadata(args: &Args) -> Result<Vec<output::InputMetadata>> {
    let mut inputs = Vec::new();
    if let Some(path) = &args.fasta_file {
        inputs.push(input_metadata("query", path)?);
    }
    if let Some(path) = &args.transcriptome_fasta {
        inputs.push(input_metadata("transcriptome", path)?);
    } else {
        inputs.push(input_metadata("dataset", &dataset::dataset_path(args)?)?);
    }
    if let Some(path) = &args.genome {
        inputs.push(input_metadata("genome", path)?);
    }
    Ok(inputs)
}

fn input_metadata(role: &str, path: &std::path::Path) -> Result<output::InputMetadata> {
    let metadata = std::fs::metadata(path)
        .with_context(|| format!("failed to read metadata for {}", path.display()))?;
    let modified_unix_seconds = metadata
        .modified()
        .ok()
        .and_then(|value| value.duration_since(UNIX_EPOCH).ok())
        .map(|value| value.as_secs());
    Ok(output::InputMetadata {
        role: role.to_string(),
        path: std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf()),
        size_bytes: metadata.len(),
        modified_unix_seconds,
    })
}

fn report_metadata(
    args: &Args,
    inputs: &[output::InputMetadata],
    phases: &[output::PhaseTiming],
    kmer_specs: &[(usize, usize)],
    total_elapsed: Duration,
) -> output::ReportMetadata {
    output::ReportMetadata {
        command: actual_command(),
        species: args.specie.clone(),
        release: args.release.clone(),
        kmer_specs: kmer_specs.to_vec(),
        hash_table_count: args.hash_table_count,
        max_on_transcriptome: args.max_on_transcriptome,
        max_on_genome: args.max_on_genome,
        threads: args.thread,
        stringent: args.stringent,
        write_kmers: args.write_kmers,
        inputs: inputs.to_vec(),
        phases: phases.to_vec(),
        total_elapsed,
        peak_rss_kib: peak_rss_kib(),
    }
}

fn actual_command() -> String {
    std::env::args_os()
        .map(|arg| shell_quote(&arg))
        .collect::<Vec<_>>()
        .join(" ")
}

fn shell_quote(value: &OsStr) -> String {
    let value = value.to_string_lossy();
    if !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"/_-.=:".contains(&byte))
    {
        value.into_owned()
    } else {
        format!("'{}'", value.replace('\'', "'\\''"))
    }
}

fn peak_rss_kib() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let value = status.lines().find(|line| line.starts_with("VmHWM:"))?;
    value.split_whitespace().nth(1)?.parse().ok()
}

fn local_dataset(args: &Args) -> Dataset {
    Dataset {
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
        transcriptome: dataset::Transcriptome::new(Vec::new()),
    }
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

fn seeds_for_occurrences<'a>(
    spec: KmerSpec,
    items: &'a [QueryItem],
    occurrences: &[Occurrence],
) -> Vec<QuerySeed<'a>> {
    occurrences
        .iter()
        .enumerate()
        .map(|(id, occ)| QuerySeed {
            id,
            seq: occurrence_seq(items, occ, spec.kmer_length),
            forward_key: occ.forward_key,
            canonical_key: occ.canonical_key,
        })
        .collect()
}

fn write_run_outputs(args: &Args, items: &[QueryItem], run: &mut SpecificRun) -> Result<()> {
    std::fs::create_dir_all(&run.output)
        .with_context(|| format!("failed to create {}", run.output.display()))?;
    write_outputs(
        args,
        items,
        &run.occurrences,
        &mut run.report,
        &run.output,
        run.spec.kmer_length,
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

fn write_outputs(
    args: &Args,
    items: &[QueryItem],
    occurrences: &[Occurrence],
    report: &mut RunReport,
    output_dir: &std::path::Path,
    kmer_length: usize,
) -> Result<()> {
    debug_assert!(
        occurrences
            .windows(2)
            .all(|pair| pair[0].item_idx <= pair[1].item_idx),
        "occurrences must be grouped by query item"
    );

    let mut kmers = args
        .write_kmers
        .then(|| output::LazyFastaWriter::new(output_dir.join("kmers.fa")));
    let mut contigs = output::LazyFastaWriter::new(output_dir.join("contigs.fa"));
    let mut masked = output::LazyFastaWriter::new(output_dir.join("masked.fa"));
    let mut occurrence_idx = 0usize;

    for (item_idx, item) in items.iter().enumerate() {
        let item_start = occurrence_idx;
        while occurrence_idx < occurrences.len() && occurrences[occurrence_idx].item_idx == item_idx
        {
            occurrence_idx += 1;
        }
        let item_occs = &occurrences[item_start..occurrence_idx];
        if item_occs.is_empty() {
            continue;
        }
        let mut specific_count = 0usize;
        let mut contig_count = 0usize;
        let mut current_contig_start = 0usize;
        let mut current_contig_no = 0usize;
        let mut last_pos = 0usize;

        for occ in item_occs {
            let seq = occurrence_seq(items, occ, kmer_length);
            if is_specific(args, item, occ) {
                specific_count += 1;
                if current_contig_start == 0 || occ.pos != last_pos + 1 {
                    if current_contig_start != 0 {
                        write_contig_record(
                            &mut contigs,
                            item,
                            current_contig_no,
                            current_contig_start,
                            last_pos,
                            kmer_length,
                        )?;
                    }
                    contig_count += 1;
                    current_contig_no = contig_count;
                    current_contig_start = occ.pos;
                }
                last_pos = occ.pos;
                if let Some(kmers) = &mut kmers {
                    write_kmer_record(kmers, item, occ, current_contig_no, seq)?;
                }
            } else {
                write_masked_record(&mut masked, item, occ, seq)?;
            }
        }

        if current_contig_start != 0 {
            write_contig_record(
                &mut contigs,
                item,
                current_contig_no,
                current_contig_start,
                last_pos,
                kmer_length,
            )?;
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

    if let Some(kmers) = &mut kmers {
        kmers.finish()?;
    }
    contigs.finish()?;
    masked.finish()?;
    Ok(())
}

fn occurrence_seq<'a>(items: &'a [QueryItem], occ: &Occurrence, kmer_length: usize) -> &'a [u8] {
    let start = occ.pos - 1;
    &items[occ.item_idx].seq[start..start + kmer_length]
}

fn contig_seq(item: &QueryItem, contig_start: usize, last_pos: usize, kmer_length: usize) -> &[u8] {
    let start = contig_start - 1;
    let end = last_pos + kmer_length - 1;
    &item.seq[start..end]
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

fn write_kmer_record(
    writer: &mut output::LazyFastaWriter,
    item: &QueryItem,
    occ: &Occurrence,
    contig_no: usize,
    seq: &[u8],
) -> Result<()> {
    writer.write_record_with(
        |writer| write_kmer_header(writer, item, occ, contig_no),
        seq,
    )
}

fn write_contig_record(
    writer: &mut output::LazyFastaWriter,
    item: &QueryItem,
    contig_no: usize,
    contig_start: usize,
    last_pos: usize,
    kmer_length: usize,
) -> Result<()> {
    let seq = contig_seq(item, contig_start, last_pos, kmer_length);
    writer.write_record_with(
        |writer| write_contig_header(writer, item, contig_no, contig_start),
        seq,
    )
}

fn write_masked_record(
    writer: &mut output::LazyFastaWriter,
    item: &QueryItem,
    occ: &Occurrence,
    seq: &[u8],
) -> Result<()> {
    writer.write_record_with(|writer| write_masked_header(writer, item, occ), seq)
}

fn write_kmer_header<W: Write>(
    writer: &mut W,
    item: &QueryItem,
    occ: &Occurrence,
    contig_no: usize,
) -> io::Result<()> {
    match item.kind {
        ItemKind::Gene => {
            write_ascii_uppercase(writer, &item.given)?;
            writer.write_all(b":")?;
            write_str(writer, item.enst.as_deref().unwrap_or("NA"))?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" ct:")?;
            write_usize(writer, contig_no)?;
            writer.write_all(b" tr:")?;
            write_u32(writer, occ.isoform_count)?;
            writer.write_all(b"/")?;
            write_u32(writer, item.isoform_total)
        }
        ItemKind::Transcript => {
            write_ascii_uppercase(writer, &item.given)?;
            writer.write_all(b":")?;
            write_str(writer, item.enst.as_deref().unwrap_or(&item.id))?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" ct:")?;
            write_usize(writer, contig_no)
        }
        ItemKind::Fasta => {
            write_str(writer, &item.id)?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" ct:")?;
            write_usize(writer, contig_no)
        }
    }
}

fn write_contig_header<W: Write>(
    writer: &mut W,
    item: &QueryItem,
    contig_no: usize,
    pos: usize,
) -> io::Result<()> {
    match item.kind {
        ItemKind::Gene => {
            write_ascii_uppercase(writer, &item.given)?;
            writer.write_all(b":")?;
            write_str(writer, item.enst.as_deref().unwrap_or("NA"))?;
            writer.write_all(b".contig_")?;
            write_usize(writer, contig_no)?;
            writer.write_all(b" (at position ")?;
            write_usize(writer, pos)?;
            writer.write_all(b")")
        }
        ItemKind::Transcript => {
            write_str(writer, item.enst.as_deref().unwrap_or(&item.id))?;
            writer.write_all(b".contig_")?;
            write_usize(writer, contig_no)?;
            writer.write_all(b" (at position ")?;
            write_usize(writer, pos)?;
            writer.write_all(b")")
        }
        ItemKind::Fasta => {
            write_str(writer, &item.id)?;
            writer.write_all(b".contig_")?;
            write_usize(writer, contig_no)?;
            writer.write_all(b" (at position ")?;
            write_usize(writer, pos)?;
            writer.write_all(b")")
        }
    }
}

fn write_masked_header<W: Write>(
    writer: &mut W,
    item: &QueryItem,
    occ: &Occurrence,
) -> io::Result<()> {
    match item.kind {
        ItemKind::Gene => {
            write_ascii_uppercase(writer, &item.given)?;
            writer.write_all(b":")?;
            write_str(writer, item.enst.as_deref().unwrap_or("NA"))?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" tr:")?;
            write_u32(writer, occ.isoform_count)?;
            writer.write_all(b"/")?;
            write_u32(writer, item.isoform_total)?;
            writer.write_all(b" genome:")?;
            write_u32(writer, occ.genome_count)?;
            writer.write_all(b" transcriptome:")?;
            write_u32(writer, occ.transcriptome_count)
        }
        ItemKind::Transcript => {
            write_ascii_uppercase(writer, &item.given)?;
            writer.write_all(b":")?;
            write_str(writer, item.enst.as_deref().unwrap_or(&item.id))?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" genome:")?;
            write_u32(writer, occ.genome_count)?;
            writer.write_all(b" transcriptome:")?;
            write_u32(writer, occ.transcriptome_count)
        }
        ItemKind::Fasta => {
            write_str(writer, &item.id)?;
            writer.write_all(b".kmer")?;
            write_usize(writer, occ.pos)?;
            writer.write_all(b" genome:")?;
            write_u32(writer, occ.genome_count)?;
            writer.write_all(b" transcriptome:")?;
            write_u32(writer, occ.transcriptome_count)
        }
    }
}

fn write_str<W: Write>(writer: &mut W, value: &str) -> io::Result<()> {
    writer.write_all(value.as_bytes())
}

fn write_usize<W: Write>(writer: &mut W, value: usize) -> io::Result<()> {
    write_u64(writer, value as u64)
}

fn write_u32<W: Write>(writer: &mut W, value: u32) -> io::Result<()> {
    write_u64(writer, value as u64)
}

fn write_u64<W: Write>(writer: &mut W, mut value: u64) -> io::Result<()> {
    let mut buffer = [0u8; 20];
    let mut idx = buffer.len();
    loop {
        idx -= 1;
        buffer[idx] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    writer.write_all(&buffer[idx..])
}

fn write_ascii_uppercase<W: Write>(writer: &mut W, value: &str) -> io::Result<()> {
    for byte in value.bytes() {
        writer.write_all(&[byte.to_ascii_uppercase()])?;
    }
    Ok(())
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
            write_kmers: true,
            thread: 2,
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
