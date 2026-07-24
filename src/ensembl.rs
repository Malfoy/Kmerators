use crate::cli::Args;
use crate::dataset::{DATASET_VERSION, Dataset, Gene, GeneInfo, TranscriptRecord, Transcriptome};
use crate::fastx;
use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::Path;

const BASE_URL: &str = "https://ftp.ensembl.org/pub";

#[derive(Debug, Clone)]
struct Meta {
    assembly: String,
    chromosomes_by_region: HashMap<String, String>,
}

#[derive(Debug, Clone)]
struct GeneRef {
    canonical_transcript_id: String,
    stable_id: String,
}

pub fn current_release(species: &str) -> Result<String> {
    let url = format!("{BASE_URL}/current/mysql/");
    let html = get_text(&url)?;
    parse_current_release(&html, species)
}

fn parse_current_release(html: &str, species: &str) -> Result<String> {
    let prefix = format!("{species}_core_");
    let link = extract_hrefs(html)
        .into_iter()
        .find(|href| href.starts_with(&prefix))
        .with_context(|| format!("species {species} not found in Ensembl current/mysql listing"))?;
    let parts = link.trim_end_matches('/').split('_').collect::<Vec<_>>();
    let release = parts
        .get(parts.len().saturating_sub(2))
        .with_context(|| format!("could not parse Ensembl release from {link}"))?;
    Ok((*release).to_string())
}

pub fn build_dataset(args: &Args) -> Result<Dataset> {
    let core_url = core_mysql_url(&args.specie, &args.release)?;
    if args.debug {
        eprintln!("Using Ensembl core URL: {core_url}");
    }

    let meta = load_meta(&core_url)?;
    let geneinfo = load_geneinfo(&core_url, &meta)?;
    let transcriptome = load_transcriptome(args, &meta)?;

    Ok(Dataset {
        species: args.specie.clone(),
        release: args.release.clone(),
        assembly: meta.assembly,
        geneinfo,
        transcriptome,
    })
}

fn core_mysql_url(species: &str, release: &str) -> Result<String> {
    let url = format!("{BASE_URL}/release-{release}/mysql/");
    let html = get_text(&url)?;
    let prefix = format!("{species}_core_{release}");
    let link = extract_hrefs(&html)
        .into_iter()
        .rfind(|href| href.starts_with(&prefix))
        .with_context(|| format!("core database for {species} release {release} not found"))?;
    Ok(format!("{url}{link}"))
}

fn load_meta(core_url: &str) -> Result<Meta> {
    let attrib_type = get_gzip_text(&format!("{core_url}attrib_type.txt.gz"))?;
    let attrib_type_id = attrib_type
        .lines()
        .filter_map(|line| {
            let fields = split_tab(line);
            (fields.get(1) == Some(&"karyotype_rank")).then(|| fields[0].to_string())
        })
        .next()
        .context("karyotype_rank attribute not found")?;

    let seq_region_attrib = get_gzip_text(&format!("{core_url}seq_region_attrib.txt.gz"))?;
    let region_ids = seq_region_attrib
        .lines()
        .filter_map(|line| {
            let fields = split_tab(line);
            (fields.get(1) == Some(&attrib_type_id.as_str())).then(|| fields[0].to_string())
        })
        .collect::<HashSet<_>>();

    let seq_region = get_gzip_text(&format!("{core_url}seq_region.txt.gz"))?;
    let mut chromosomes_by_region = HashMap::new();
    for line in seq_region.lines() {
        let fields = split_tab(line);
        if fields.len() >= 2 && region_ids.contains(fields[0]) {
            chromosomes_by_region.insert(fields[0].to_string(), fields[1].to_string());
        }
    }

    let meta = get_gzip_text(&format!("{core_url}meta.txt.gz"))?;
    let assembly = meta
        .lines()
        .filter_map(|line| {
            let fields = split_tab(line);
            (fields.get(2) == Some(&"assembly.default")).then(|| {
                fields
                    .get(3)
                    .copied()
                    .unwrap_or("unknown")
                    .split('.')
                    .next()
                    .unwrap_or("unknown")
                    .to_string()
            })
        })
        .next()
        .context("assembly.default not found in Ensembl meta")?;

    Ok(Meta {
        assembly,
        chromosomes_by_region,
    })
}

fn load_geneinfo(core_url: &str, meta: &Meta) -> Result<GeneInfo> {
    let gene_text = get_gzip_text(&format!("{core_url}gene.txt.gz"))?;
    let mut genes = HashMap::<String, Gene>::new();
    let mut gene_refs = HashMap::<String, GeneRef>::new();
    let mut xrefs = HashMap::<String, Vec<String>>::new();

    for line in gene_text.lines() {
        let fields = split_tab(line);
        if fields.len() <= 12 {
            continue;
        }
        let Some(chr) = meta.chromosomes_by_region.get(fields[3]) else {
            continue;
        };
        let stable_id = fields[12].to_string();
        gene_refs.insert(
            fields[0].to_string(),
            GeneRef {
                canonical_transcript_id: fields[11].to_string(),
                stable_id: stable_id.clone(),
            },
        );
        xrefs
            .entry(fields[7].to_string())
            .or_default()
            .push(fields[0].to_string());
        genes.insert(
            stable_id,
            Gene {
                symbol: None,
                aliases: Vec::new(),
                canonical: None,
                transcripts: Vec::new(),
                chr: chr.clone(),
                start: fields[4].parse().unwrap_or(0),
                end: fields[5].parse().unwrap_or(0),
                strand: fields[6].parse().unwrap_or(0),
                biotype: fields[1].to_string(),
                desc: fields[9].to_string(),
            },
        );
    }

    let xref_text = get_gzip_text(&format!("{core_url}xref.txt.gz"))?;
    let mut symbols = HashMap::<String, Vec<String>>::new();
    for line in xref_text.lines() {
        let fields = split_tab(line);
        if fields.len() <= 3 {
            continue;
        }
        let Some(gene_ids) = xrefs.get(fields[0]) else {
            continue;
        };
        for gene_id in gene_ids {
            let Some(gene_ref) = gene_refs.get(gene_id) else {
                continue;
            };
            if let Some(gene) = genes.get_mut(&gene_ref.stable_id) {
                let symbol = fields[3].to_string();
                if symbol != "\\N" {
                    gene.symbol = Some(symbol.clone());
                    symbols
                        .entry(symbol.to_ascii_uppercase())
                        .or_default()
                        .push(gene_ref.stable_id.clone());
                }
            }
        }
    }

    let synonym_text = get_gzip_text(&format!("{core_url}external_synonym.txt.gz"))?;
    let mut alias_sets = HashMap::<String, HashSet<String>>::new();
    for line in synonym_text.lines() {
        let fields = split_tab(line);
        if fields.len() < 2 {
            continue;
        }
        let Some(gene_ids) = xrefs.get(fields[0]) else {
            continue;
        };
        for gene_id in gene_ids {
            let Some(gene_ref) = gene_refs.get(gene_id) else {
                continue;
            };
            if let Some(gene) = genes.get_mut(&gene_ref.stable_id) {
                let alias = fields[1].to_string();
                gene.aliases.push(alias.clone());
                alias_sets
                    .entry(alias.to_ascii_uppercase())
                    .or_default()
                    .insert(gene_ref.stable_id.clone());
            }
        }
    }
    let aliases = alias_sets
        .into_iter()
        .map(|(key, set)| {
            let mut vals = set.into_iter().collect::<Vec<_>>();
            vals.sort();
            (key, vals)
        })
        .collect::<HashMap<_, _>>();

    let transcript_text = get_gzip_text(&format!("{core_url}transcript.txt.gz"))?;
    let mut transcripts = HashMap::<String, String>::new();
    for line in transcript_text.lines() {
        let fields = split_tab(line);
        if fields.len() <= 13 {
            continue;
        }
        let Some(gene_ref) = gene_refs.get(fields[1]) else {
            continue;
        };
        let enst = fields[13].to_string();
        transcripts.insert(enst.to_ascii_uppercase(), gene_ref.stable_id.clone());
        if let Some(gene) = genes.get_mut(&gene_ref.stable_id) {
            gene.transcripts.push(enst.clone());
            if fields[0] == gene_ref.canonical_transcript_id {
                gene.canonical = Some(enst);
            }
        }
    }

    Ok(GeneInfo {
        assembly: meta.assembly.clone(),
        chromosomes: meta.chromosomes_by_region.values().cloned().collect(),
        version: DATASET_VERSION,
        genes,
        symbols,
        aliases,
        transcripts,
    })
}

fn load_transcriptome(args: &Args, meta: &Meta) -> Result<Transcriptome> {
    let chromosomes = meta
        .chromosomes_by_region
        .values()
        .cloned()
        .collect::<HashSet<_>>();
    let mut records_by_id = HashMap::<String, Vec<u8>>::new();

    for kind in ["cdna", "ncrna"] {
        let dir_url = format!(
            "{BASE_URL}/release-{}/fasta/{}/{kind}/",
            args.release, args.specie
        );
        let html = get_text(&dir_url)?;
        let pattern = if kind == "cdna" {
            "cdna.all.fa.gz"
        } else {
            "ncrna.fa.gz"
        };
        let file = extract_hrefs(&html)
            .into_iter()
            .find(|href| href.contains(pattern))
            .with_context(|| format!("could not find {pattern} in {dir_url}"))?;
        let url = format!("{dir_url}{file}");
        if args.debug {
            eprintln!("Downloading transcriptome FASTA: {url}");
        }
        load_transcriptome_fasta(&url, &chromosomes, &mut records_by_id)?;
    }

    let mut records = records_by_id
        .into_iter()
        .map(|(id, seq)| TranscriptRecord { id, seq })
        .collect::<Vec<_>>();
    records.sort_by(|a, b| a.id.cmp(&b.id));
    Ok(Transcriptome::new(records))
}

fn load_transcriptome_fasta(
    url: &str,
    chromosomes: &HashSet<String>,
    records_by_id: &mut HashMap<String, Vec<u8>>,
) -> Result<()> {
    let bytes = get_bytes(url)?;
    let tmp = tempfile::Builder::new()
        .suffix(".fa.gz")
        .tempfile()
        .context("failed to create temporary compressed transcriptome FASTA")?;
    std::fs::write(tmp.path(), &bytes).with_context(|| {
        format!(
            "failed to write temporary compressed FASTA {}",
            tmp.path().display()
        )
    })?;

    for rec in fastx::read_records(tmp.path())? {
        let fields = rec
            .header
            .split(|b| b.is_ascii_whitespace())
            .collect::<Vec<_>>();
        if fields.len() <= 2 {
            continue;
        }
        let location = String::from_utf8_lossy(fields[2]);
        let chr = location.split(':').nth(2).unwrap_or("");
        if !chromosomes.contains(chr) {
            continue;
        }
        let id = fastx::first_token_without_version(&rec.header);
        records_by_id.entry(id).or_insert(rec.seq);
    }
    Ok(())
}

fn get_text(url: &str) -> Result<String> {
    let response = reqwest::blocking::get(url).with_context(|| format!("GET {url} failed"))?;
    if !response.status().is_success() {
        bail!("GET {url} returned {}", response.status());
    }
    response
        .text()
        .with_context(|| format!("failed to read response body from {url}"))
}

fn get_bytes(url: &str) -> Result<Vec<u8>> {
    let response = reqwest::blocking::get(url).with_context(|| format!("GET {url} failed"))?;
    if !response.status().is_success() {
        bail!("GET {url} returned {}", response.status());
    }
    response
        .bytes()
        .map(|b| b.to_vec())
        .with_context(|| format!("failed to read response body from {url}"))
}

fn get_gzip_text(url: &str) -> Result<String> {
    let bytes = get_bytes(url)?;
    let mut decoder = GzDecoder::new(bytes.as_slice());
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .with_context(|| format!("failed to decompress {url}"))?;
    Ok(text)
}

fn split_tab(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn extract_hrefs(html: &str) -> Vec<String> {
    let mut hrefs = Vec::new();
    let mut rest = html;
    while let Some(pos) = rest.find("href=\"") {
        rest = &rest[pos + 6..];
        let Some(end) = rest.find('"') else {
            break;
        };
        let href = &rest[..end];
        if !href.starts_with('?') && href != "../" {
            hrefs.push(href.to_string());
        }
        rest = &rest[end + 1..];
    }
    hrefs
}

#[allow(dead_code)]
fn path_exists(path: &Path) -> bool {
    path.exists()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_species_release_from_current_mysql_listing() {
        let html = r#"
            <a href="homo_sapiens_core_116_38/">homo_sapiens_core_116_38/</a>
            <a href="mus_musculus_core_116_39/">mus_musculus_core_116_39/</a>
        "#;

        assert_eq!(parse_current_release(html, "homo_sapiens").unwrap(), "116");
        assert_eq!(parse_current_release(html, "mus_musculus").unwrap(), "116");
    }

    #[test]
    fn rejects_species_missing_from_current_mysql_listing() {
        let err = parse_current_release(
            r#"<a href="homo_sapiens_core_116_38/">human</a>"#,
            "toy_species",
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("species toy_species not found in Ensembl current/mysql listing")
        );
    }
}
