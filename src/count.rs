use crate::fastx;
use crate::kmer::{self, KeyMode};
use anyhow::{Context, Result};
use crossbeam_channel::bounded;
use hashbrown::HashMap;
use helicase::input::*;
use helicase::*;
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::{PackedSeqVec, Seq, SeqVec};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

const FASTX_CONFIG: Config = ParserOptions::default().config();
const BATCH_BASES: usize = 8 * 1024 * 1024;
const SCAN_BLOCK_BASES: usize = 4 * 1024 * 1024;

type KmerId = usize;
type Partition = HashMap<u64, Vec<KmerId>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountKey {
    Forward,
    Canonical,
}

#[derive(Debug)]
pub struct QuerySeed<'a> {
    pub id: usize,
    pub seq: &'a [u8],
    pub forward_key: u64,
    pub canonical_key: u64,
}

#[derive(Debug)]
pub struct CountIndex {
    k: usize,
    m: usize,
    key_mode: KeyMode,
    key_kind: CountKey,
    partitions: Vec<Partition>,
    active_simd_partitions: Vec<bool>,
    query_count: usize,
    counter_members: Vec<Vec<KmerId>>,
}

impl CountIndex {
    pub fn new(
        k: usize,
        m: usize,
        hash_table_count: usize,
        key_kind: CountKey,
        seeds: &[QuerySeed<'_>],
    ) -> Self {
        let key_mode = KeyMode::for_k(k);
        let mut index = Self {
            k,
            m,
            key_mode,
            key_kind,
            partitions: (0..hash_table_count).map(|_| Partition::new()).collect(),
            active_simd_partitions: vec![false; hash_table_count],
            query_count: seeds.len(),
            counter_members: Vec::new(),
        };
        let mut counter_by_key = HashMap::<u64, KmerId>::new();

        for seed in seeds {
            let key = match key_kind {
                CountKey::Forward => seed.forward_key,
                CountKey::Canonical => seed.canonical_key,
            };
            let counter_id = if let Some(&counter_id) = counter_by_key.get(&key) {
                index.counter_members[counter_id].push(seed.id);
                counter_id
            } else {
                let counter_id = index.counter_members.len();
                counter_by_key.insert(key, counter_id);
                index.counter_members.push(vec![seed.id]);
                counter_id
            };

            match key_kind {
                CountKey::Forward => {
                    if let Some(minimizer) = scalar_minimizer(seed.seq, m) {
                        index.insert(minimizer, seed.forward_key, counter_id);
                    }
                    index.insert_simd_minimizers(seed.seq);
                }
                CountKey::Canonical => {
                    if let Some(minimizer) = scalar_minimizer(seed.seq, m) {
                        index.insert(minimizer, seed.canonical_key, counter_id);
                    }
                    index.insert_simd_minimizers(seed.seq);
                    if let Some(rc) = kmer::revcomp(seed.seq) {
                        if let Some(minimizer) = scalar_minimizer(&rc, m) {
                            index.insert(minimizer, seed.canonical_key, counter_id);
                        }
                        index.insert_simd_minimizers(&rc);
                    }
                }
            }
        }

        index
    }

    fn insert(&mut self, minimizer: u64, key: u64, id: KmerId) {
        let table_id = table_for_minimizer(minimizer, self.partitions.len());
        let ids = self.partitions[table_id].entry(key).or_default();
        if !ids.contains(&id) {
            ids.push(id);
        }
    }

    fn insert_simd_minimizers(&mut self, seq: &[u8]) {
        if seq.len() < self.k {
            return;
        }
        let window = self.k - self.m + 1;
        let mut packed = PackedSeqVec::default();
        packed.push_ascii(seq);
        let builder = minimizers(self.m, window);
        let mut positions = Vec::new();
        builder.run(packed.as_slice(), &mut positions);
        for pos in positions {
            let pos = pos as usize;
            let table_id = table_for_minimizer(
                packed.slice(pos..pos + self.m).as_u64(),
                self.active_simd_partitions.len(),
            );
            self.active_simd_partitions[table_id] = true;
        }
    }

    fn counter_count(&self) -> usize {
        self.counter_members.len()
    }

    fn expand_counts(&self, compact: Vec<u32>) -> Vec<u32> {
        debug_assert_eq!(compact.len(), self.counter_members.len());
        let mut expanded = vec![0u32; self.query_count];
        for (count, members) in compact.into_iter().zip(&self.counter_members) {
            for &query_id in members {
                expanded[query_id] = count;
            }
        }
        expanded
    }

    pub fn is_empty(&self) -> bool {
        self.partitions.iter().all(HashMap::is_empty)
    }
}

#[allow(dead_code)]
pub fn count_path(path: &Path, index: Arc<CountIndex>, threads: usize) -> Result<Vec<u32>> {
    Ok(count_path_many(path, vec![index], threads)?.remove(0))
}

pub fn count_path_many(
    path: &Path,
    indexes: Vec<Arc<CountIndex>>,
    threads: usize,
) -> Result<Vec<Vec<u32>>> {
    let shapes = index_shapes(&indexes);
    let active = active_indexes(&indexes);
    if active.is_empty() {
        return Ok(expand_counts(&indexes, zero_counts(&shapes)));
    }

    let active = Arc::new(active);
    let (tx, rx) = bounded::<Vec<Vec<u8>>>(threads * 2);
    let path = PathBuf::from(path);
    let producer = thread::spawn(move || produce_batches(&path, tx));

    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = rx.clone();
        let active = Arc::clone(&active);
        let shapes = shapes.clone();
        handles.push(thread::spawn(move || {
            let mut counts = zero_counts(&shapes);
            while let Ok(batch) = rx.recv() {
                for seq in batch {
                    for &(idx, ref index) in active.iter() {
                        scan_sequence(&seq, index, &mut counts[idx]);
                    }
                }
            }
            counts
        }));
    }

    producer
        .join()
        .expect("producer thread panicked")
        .context("failed while reading reference sequences")?;

    let mut total = zero_counts(&shapes);
    for handle in handles {
        let local = handle.join().expect("worker thread panicked");
        add_counts(&mut total, local);
    }

    Ok(expand_counts(&indexes, total))
}

#[allow(dead_code)]
pub fn count_transcriptome_records(
    records: &[crate::dataset::TranscriptRecord],
    index: Arc<CountIndex>,
    threads: usize,
) -> Vec<u32> {
    count_transcriptome_records_many(records, vec![index], threads).remove(0)
}

pub fn count_transcriptome_records_many(
    records: &[crate::dataset::TranscriptRecord],
    indexes: Vec<Arc<CountIndex>>,
    threads: usize,
) -> Vec<Vec<u32>> {
    let shapes = index_shapes(&indexes);
    let active = active_indexes(&indexes);
    if active.is_empty() {
        return expand_counts(&indexes, zero_counts(&shapes));
    }

    let active = Arc::new(active);

    let chunk_size = (records.len() / threads.max(1)).max(1);
    let mut handles = Vec::new();
    for chunk in records.chunks(chunk_size) {
        let seqs = chunk.iter().map(|rec| rec.seq.clone()).collect::<Vec<_>>();
        let active = Arc::clone(&active);
        let shapes = shapes.clone();
        handles.push(thread::spawn(move || {
            let mut counts = zero_counts(&shapes);
            for seq in seqs {
                for &(idx, ref index) in active.iter() {
                    scan_sequence(&seq, index, &mut counts[idx]);
                }
            }
            counts
        }));
    }

    let mut total = zero_counts(&shapes);
    for handle in handles {
        let local = handle.join().expect("transcriptome worker panicked");
        add_counts(&mut total, local);
    }
    expand_counts(&indexes, total)
}

fn index_shapes(indexes: &[Arc<CountIndex>]) -> Vec<usize> {
    indexes.iter().map(|index| index.counter_count()).collect()
}

fn expand_counts(indexes: &[Arc<CountIndex>], compact: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    indexes
        .iter()
        .zip(compact)
        .map(|(index, counts)| index.expand_counts(counts))
        .collect()
}

fn active_indexes(indexes: &[Arc<CountIndex>]) -> Vec<(usize, Arc<CountIndex>)> {
    indexes
        .iter()
        .enumerate()
        .filter(|(_, index)| !index.is_empty())
        .map(|(idx, index)| (idx, Arc::clone(index)))
        .collect()
}

fn zero_counts(shapes: &[usize]) -> Vec<Vec<u32>> {
    shapes.iter().map(|&count| vec![0u32; count]).collect()
}

fn add_counts(total: &mut [Vec<u32>], local: Vec<Vec<u32>>) {
    for (total_run, local_run) in total.iter_mut().zip(local) {
        for (dst, src) in total_run.iter_mut().zip(local_run) {
            *dst = dst.saturating_add(src);
        }
    }
}

fn produce_batches(path: &Path, tx: crossbeam_channel::Sender<Vec<Vec<u8>>>) -> Result<()> {
    let mut parser = FastxParser::<FASTX_CONFIG>::from_file(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
    let mut batch = Vec::new();
    let mut bases = 0usize;

    while let Some(_) = parser.next() {
        let seq = fastx::normalize_dna(parser.get_dna_string());
        if seq.is_empty() {
            continue;
        }
        bases += seq.len();
        batch.push(seq);
        if bases >= BATCH_BASES {
            if tx.send(std::mem::take(&mut batch)).is_err() {
                return Ok(());
            }
            bases = 0;
        }
    }

    if !batch.is_empty() {
        let _ = tx.send(batch);
    }
    Ok(())
}

fn scan_sequence(seq: &[u8], index: &CountIndex, counts: &mut [u32]) {
    for (_, chunk) in kmer::valid_chunks(seq) {
        if chunk.len() < index.k {
            continue;
        }
        let mut start = 0usize;
        while start + index.k <= chunk.len() {
            let remaining_windows = chunk.len() - start - index.k + 1;
            let window_limit = remaining_windows.min(SCAN_BLOCK_BASES);
            let end = (start + window_limit + index.k - 1).min(chunk.len());
            let block = &chunk[start..end];
            if has_any_simd_minimizer(block, index) {
                scan_valid_chunk(block, index, counts, window_limit);
            }
            start += window_limit;
        }
    }
}

fn scan_valid_chunk(seq: &[u8], index: &CountIndex, counts: &mut [u32], window_limit: usize) {
    let kmer_count = seq.len() - index.k + 1;
    let kmer_count = kmer_count.min(window_limit);
    let m_values = kmer::minimizer_values(seq, index.m);
    let minimizers = kmer::window_minima(&m_values, index.k - index.m + 1);
    debug_assert!(minimizers.len() >= kmer_count);

    match index.key_mode {
        KeyMode::ExactU64 => scan_exact_chunk(seq, index, counts, &minimizers[..kmer_count]),
        KeyMode::HashedU64 => {
            for pos in 0..kmer_count {
                let minimizer = minimizers[pos];
                let partition =
                    &index.partitions[table_for_minimizer(minimizer, index.partitions.len())];
                if partition.is_empty() {
                    continue;
                }
                let kmer_seq = &seq[pos..pos + index.k];
                let key = match index.key_kind {
                    CountKey::Forward => kmer::forward_key(kmer_seq, index.key_mode),
                    CountKey::Canonical => kmer::canonical_key(kmer_seq, index.key_mode),
                };
                if let Some(ids) = key.and_then(|key| partition.get(&key)) {
                    for &id in ids {
                        counts[id] = counts[id].saturating_add(1);
                    }
                }
            }
        }
    }
}

fn scan_exact_chunk(seq: &[u8], index: &CountIndex, counts: &mut [u32], minimizers: &[u64]) {
    let keys = kmer::kmer_keys_for_chunk(seq, index.k, KeyMode::ExactU64);
    for (pos, (forward, canonical)) in keys.into_iter().take(minimizers.len()).enumerate() {
        let minimizer = minimizers[pos];
        let partition = &index.partitions[table_for_minimizer(minimizer, index.partitions.len())];
        if partition.is_empty() {
            continue;
        }
        let key = match index.key_kind {
            CountKey::Forward => forward,
            CountKey::Canonical => canonical,
        };
        if let Some(ids) = partition.get(&key) {
            for &id in ids {
                counts[id] = counts[id].saturating_add(1);
            }
        }
    }
}

fn scalar_minimizer(seq: &[u8], m: usize) -> Option<u64> {
    if seq.len() < m {
        return None;
    }
    kmer::minimizer_values(seq, m).into_iter().min()
}

fn table_for_minimizer(minimizer: u64, table_count: usize) -> usize {
    if table_count.is_power_of_two() {
        (minimizer as usize) & (table_count - 1)
    } else {
        (minimizer as usize) % table_count
    }
}

fn has_any_simd_minimizer(seq: &[u8], index: &CountIndex) -> bool {
    if seq.len() < index.k || !index.active_simd_partitions.iter().any(|&active| active) {
        return false;
    }

    let window = index.k - index.m + 1;
    let mut packed = PackedSeqVec::default();
    packed.push_ascii(seq);
    let builder = minimizers(index.m, window);
    let mut positions = Vec::new();
    builder.run(packed.as_slice(), &mut positions);

    positions.iter().copied().any(|pos| {
        let pos = pos as usize;
        let value = packed.slice(pos..pos + index.m).as_u64();
        index.active_simd_partitions[table_for_minimizer(value, index.active_simd_partitions.len())]
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::TranscriptRecord;
    use std::sync::Arc;

    fn seed<'a>(id: usize, seq: &'a [u8]) -> QuerySeed<'a> {
        let mode = KeyMode::for_k(seq.len());
        QuerySeed {
            id,
            seq,
            forward_key: kmer::forward_key(seq, mode).unwrap(),
            canonical_key: kmer::canonical_key(seq, mode).unwrap(),
        }
    }

    #[test]
    fn canonical_counts_include_reverse_complements() {
        let records = vec![TranscriptRecord {
            id: "TR1".to_string(),
            seq: b"CGTACG".to_vec(),
        }];
        let seeds = vec![seed(0, b"ACG")];

        let forward = Arc::new(CountIndex::new(3, 2, 7, CountKey::Forward, &seeds));
        let canonical = Arc::new(CountIndex::new(3, 2, 7, CountKey::Canonical, &seeds));

        assert_eq!(count_transcriptome_records(&records, forward, 1), vec![1]);
        assert_eq!(count_transcriptome_records(&records, canonical, 1), vec![2]);
    }

    #[test]
    fn ambiguous_bases_split_reference_chunks() {
        let records = vec![TranscriptRecord {
            id: "TR1".to_string(),
            seq: b"ACGNACG".to_vec(),
        }];
        let seeds = vec![seed(0, b"ACG")];
        let index = Arc::new(CountIndex::new(3, 2, 8, CountKey::Forward, &seeds));

        assert_eq!(count_transcriptome_records(&records, index, 2), vec![2]);
    }

    #[test]
    fn multi_count_matches_single_indexes() {
        let records = vec![TranscriptRecord {
            id: "TR1".to_string(),
            seq: b"ACGTACCCC".to_vec(),
        }];
        let seeds_k3 = vec![seed(0, b"ACG"), seed(1, b"CCC")];
        let seeds_k4 = vec![seed(0, b"ACGT")];

        let single_k3 = Arc::new(CountIndex::new(3, 2, 8, CountKey::Forward, &seeds_k3));
        let single_k4 = Arc::new(CountIndex::new(4, 2, 8, CountKey::Forward, &seeds_k4));
        let many_k3 = Arc::new(CountIndex::new(3, 2, 8, CountKey::Forward, &seeds_k3));
        let many_k4 = Arc::new(CountIndex::new(4, 2, 8, CountKey::Forward, &seeds_k4));

        let expected_k3 = count_transcriptome_records(&records, single_k3, 2);
        let expected_k4 = count_transcriptome_records(&records, single_k4, 2);
        let actual = count_transcriptome_records_many(&records, vec![many_k3, many_k4], 2);

        assert_eq!(actual, vec![expected_k3, expected_k4]);
    }

    #[test]
    fn duplicate_query_kmers_share_a_counter_and_expand_to_each_position() {
        let records = vec![TranscriptRecord {
            id: "TR1".to_string(),
            seq: b"ACGACG".to_vec(),
        }];
        let seeds = vec![seed(0, b"ACG"), seed(1, b"ACG")];
        let index = Arc::new(CountIndex::new(3, 2, 8, CountKey::Forward, &seeds));

        assert_eq!(index.counter_count(), 1);
        assert_eq!(count_transcriptome_records(&records, index, 2), vec![2, 2]);
    }
}
