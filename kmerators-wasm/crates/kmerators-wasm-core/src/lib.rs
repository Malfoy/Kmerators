use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::VecDeque;

use rustc_hash::FxHashMap as HashMap;

const MAX_EXACT_K: usize = 31;
const MAX_MINIMIZER_K: usize = 31;
const PREFIX_FILTER_MAX_BITS: usize = 20;
const HASH_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const HASH_PRIME: u64 = 0x0000_0100_0000_01b3;

thread_local! {
    static SESSIONS: RefCell<Vec<Option<KmeratorRun>>> = const { RefCell::new(Vec::new()) };
    static LAST_ERROR: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[derive(Debug, Clone, Copy)]
pub struct RunParams {
    pub k: usize,
    pub minimizer_length: usize,
    pub max_transcriptome_count: u32,
    pub max_genome_count: u32,
}

#[derive(Debug)]
pub struct KmeratorRun {
    params: RunParams,
    query_items: Vec<QueryItem>,
    occurrences: Vec<Occurrence>,
    transcriptome_index: QueryIndex,
    genome_index: QueryIndex,
    transcriptome_counts: Vec<u32>,
    genome_counts: Vec<u32>,
    transcriptome_enabled: bool,
    genome_enabled: bool,
    transcriptome_scanner: FastxStreamScanner,
    genome_scanner: FastxStreamScanner,
    last_json: Vec<u8>,
}

#[derive(Debug, Clone)]
struct QueryItem {
    id: String,
    seq: Vec<u8>,
}

#[derive(Debug, Clone)]
struct Occurrence {
    item_idx: usize,
    pos: usize,
    forward_key: u64,
    canonical_key: u64,
    minimizer: u64,
    revcomp_minimizer: u64,
}

#[derive(Debug, Clone)]
struct QueryIndex {
    map: HashMap<u64, HashMap<u64, Vec<usize>>>,
    filter: PrefixFilter,
}

#[derive(Debug, Clone)]
struct PrefixFilter {
    k: usize,
    bits: usize,
    words: Vec<u64>,
}

#[derive(Debug, Clone, Copy)]
enum KeyMode {
    Exact,
    Hashed,
}

#[derive(Debug, Clone, Copy)]
enum CountKey {
    Forward,
    Canonical,
}

#[derive(Debug)]
struct FastxStreamScanner {
    mode: StreamMode,
    line: Vec<u8>,
    fastq_state: FastqState,
    fastq_seq_len: usize,
    fastq_quality_len: usize,
    roller: KmerRoller,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StreamMode {
    Unknown,
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastqState {
    Header,
    Sequence,
    Quality,
}

#[derive(Debug)]
struct KmerRoller {
    k: usize,
    m: usize,
    mode: KeyMode,
    valid_len: usize,
    forward: u64,
    reverse: u64,
    exact_mask: u64,
    rc_shift: usize,
    minimizer_forward: u64,
    minimizer_mask: u64,
    minimizers: VecDeque<(usize, u64)>,
    window: VecDeque<u8>,
}

#[derive(Debug, Clone, Copy)]
struct KmerHit {
    forward_key: u64,
    canonical_key: u64,
    minimizer: u64,
}

#[derive(Debug, Deserialize, Serialize)]
struct RunOutput {
    kmer_length: usize,
    minimizer_length: usize,
    query_sequences: usize,
    query_kmers: usize,
    retained_kmers: usize,
    masked_kmers: usize,
    transcriptome_enabled: bool,
    genome_enabled: bool,
    transcriptome_hits: Option<u64>,
    genome_hits: Option<u64>,
    files: OutputFiles,
}

#[derive(Debug, Deserialize, Serialize)]
struct OutputFiles {
    contigs_fa: String,
    masked_fa: String,
    report_md: String,
}

impl RunParams {
    pub fn validate(self) -> Result<Self, String> {
        if self.k == 0 {
            return Err("k-mer length must be greater than zero".to_string());
        }
        let minimizer_length = if self.minimizer_length == 0 {
            9.min(self.k)
        } else {
            self.minimizer_length
        };
        if minimizer_length == 0 || minimizer_length > self.k {
            return Err(format!(
                "minimizer length must be in 1..={} for k={}",
                self.k, self.k
            ));
        }
        if minimizer_length > MAX_MINIMIZER_K {
            return Err(format!(
                "minimizer length must be <= {MAX_MINIMIZER_K} in the browser build"
            ));
        }
        Ok(Self {
            minimizer_length,
            ..self
        })
    }
}

impl KmeratorRun {
    pub fn new(params: RunParams) -> Result<Self, String> {
        let params = params.validate()?;
        Ok(Self {
            params,
            query_items: Vec::new(),
            occurrences: Vec::new(),
            transcriptome_index: QueryIndex::new(params.k),
            genome_index: QueryIndex::new(params.k),
            transcriptome_counts: Vec::new(),
            genome_counts: Vec::new(),
            transcriptome_enabled: false,
            genome_enabled: false,
            transcriptome_scanner: FastxStreamScanner::new(params.k, params.minimizer_length),
            genome_scanner: FastxStreamScanner::new(params.k, params.minimizer_length),
            last_json: Vec::new(),
        })
    }

    pub fn set_query_fasta(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.query_items = parse_query_fastx(bytes)?;
        if self.query_items.is_empty() {
            return Err("query FASTA/FASTQ contains no usable sequence".to_string());
        }

        self.occurrences = build_occurrences(
            self.params.k,
            self.params.minimizer_length,
            &self.query_items,
        );
        if self.occurrences.is_empty() {
            return Err("query FASTA/FASTQ produced no valid k-mers".to_string());
        }

        self.transcriptome_index = QueryIndex::with_capacity(self.params.k, self.occurrences.len());
        self.genome_index = QueryIndex::with_capacity(self.params.k, self.occurrences.len());
        for (idx, occ) in self.occurrences.iter().enumerate() {
            self.transcriptome_index
                .insert(occ.minimizer, occ.forward_key, idx);
            self.genome_index
                .insert(occ.minimizer, occ.canonical_key, idx);
            if occ.revcomp_minimizer != occ.minimizer {
                self.genome_index
                    .insert(occ.revcomp_minimizer, occ.canonical_key, idx);
            }
        }

        self.transcriptome_counts = vec![0; self.occurrences.len()];
        self.genome_counts = vec![0; self.occurrences.len()];
        self.transcriptome_scanner =
            FastxStreamScanner::new(self.params.k, self.params.minimizer_length);
        self.genome_scanner = FastxStreamScanner::new(self.params.k, self.params.minimizer_length);
        self.last_json.clear();
        Ok(())
    }

    pub fn set_transcriptome_enabled(&mut self, enabled: bool) {
        self.transcriptome_enabled = enabled;
    }

    pub fn set_genome_enabled(&mut self, enabled: bool) {
        self.genome_enabled = enabled;
    }

    pub fn scan_transcriptome_chunk(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.ensure_query_ready()?;
        self.transcriptome_enabled = true;
        scan_chunk(
            bytes,
            CountKey::Forward,
            &self.transcriptome_index,
            &mut self.transcriptome_counts,
            &mut self.transcriptome_scanner,
        );
        Ok(())
    }

    pub fn finish_transcriptome_source(&mut self) {
        self.transcriptome_scanner.finish_source(
            CountKey::Forward,
            &self.transcriptome_index,
            &mut self.transcriptome_counts,
        );
    }

    pub fn scan_genome_chunk(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.ensure_query_ready()?;
        self.genome_enabled = true;
        scan_chunk(
            bytes,
            CountKey::Canonical,
            &self.genome_index,
            &mut self.genome_counts,
            &mut self.genome_scanner,
        );
        Ok(())
    }

    pub fn finish_genome_source(&mut self) {
        self.genome_scanner.finish_source(
            CountKey::Canonical,
            &self.genome_index,
            &mut self.genome_counts,
        );
    }

    pub fn finish(&mut self) -> Result<&[u8], String> {
        self.ensure_query_ready()?;
        let output = self.build_output();
        self.last_json =
            serde_json::to_vec(&output).map_err(|err| format!("failed to encode output: {err}"))?;
        Ok(&self.last_json)
    }

    pub fn last_json(&self) -> &[u8] {
        &self.last_json
    }

    fn ensure_query_ready(&self) -> Result<(), String> {
        if self.occurrences.is_empty() {
            Err("load a query FASTA/FASTQ before scanning references".to_string())
        } else {
            Ok(())
        }
    }

    fn build_output(&self) -> RunOutput {
        let retained = self.retained_flags();
        let retained_kmers = retained.iter().filter(|&&keep| keep).count();
        let masked_kmers = retained.len() - retained_kmers;
        let files = OutputFiles {
            contigs_fa: self.write_contigs_fasta(&retained),
            masked_fa: self.write_masked_fasta(&retained),
            report_md: self.write_report(retained_kmers, masked_kmers),
        };
        RunOutput {
            kmer_length: self.params.k,
            minimizer_length: self.params.minimizer_length,
            query_sequences: self.query_items.len(),
            query_kmers: self.occurrences.len(),
            retained_kmers,
            masked_kmers,
            transcriptome_enabled: self.transcriptome_enabled,
            genome_enabled: self.genome_enabled,
            transcriptome_hits: self.transcriptome_enabled.then(|| {
                self.transcriptome_counts
                    .iter()
                    .map(|&count| u64::from(count))
                    .sum()
            }),
            genome_hits: self.genome_enabled.then(|| {
                self.genome_counts
                    .iter()
                    .map(|&count| u64::from(count))
                    .sum()
            }),
            files,
        }
    }

    fn retained_flags(&self) -> Vec<bool> {
        self.occurrences
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                (!self.transcriptome_enabled
                    || self.transcriptome_counts[idx] <= self.params.max_transcriptome_count)
                    && (!self.genome_enabled
                        || self.genome_counts[idx] <= self.params.max_genome_count)
            })
            .collect()
    }
    fn write_masked_fasta(&self, retained: &[bool]) -> String {
        let mut text = String::new();
        for (idx, occ) in self.occurrences.iter().enumerate() {
            if retained[idx] {
                continue;
            }
            let item = &self.query_items[occ.item_idx];
            let seq = &item.seq[occ.pos..occ.pos + self.params.k];
            push_fasta_record(
                &mut text,
                &format!(
                    "{}:{}-{} tr={} ge={}",
                    item.id,
                    occ.pos + 1,
                    occ.pos + self.params.k,
                    count_label(self.transcriptome_enabled, self.transcriptome_counts[idx]),
                    count_label(self.genome_enabled, self.genome_counts[idx])
                ),
                seq,
            );
        }
        text
    }

    fn write_contigs_fasta(&self, retained: &[bool]) -> String {
        let mut text = String::new();
        let mut current_item = None::<usize>;
        let mut start = 0usize;
        let mut last = 0usize;
        let mut contig_idx = 1usize;

        for (occ_idx, occ) in self.occurrences.iter().enumerate() {
            if !retained[occ_idx] {
                continue;
            }

            match current_item {
                Some(item_idx) if item_idx == occ.item_idx && occ.pos == last + 1 => {
                    last = occ.pos;
                }
                Some(item_idx) if item_idx == occ.item_idx => {
                    push_contig(
                        &mut text,
                        &self.query_items[item_idx],
                        contig_idx,
                        start,
                        last,
                        self.params.k,
                    );
                    contig_idx += 1;
                    start = occ.pos;
                    last = occ.pos;
                }
                Some(item_idx) => {
                    push_contig(
                        &mut text,
                        &self.query_items[item_idx],
                        contig_idx,
                        start,
                        last,
                        self.params.k,
                    );
                    current_item = Some(occ.item_idx);
                    contig_idx = 1;
                    start = occ.pos;
                    last = occ.pos;
                }
                None => {
                    current_item = Some(occ.item_idx);
                    start = occ.pos;
                    last = occ.pos;
                }
            }
        }

        if let Some(item_idx) = current_item {
            push_contig(
                &mut text,
                &self.query_items[item_idx],
                contig_idx,
                start,
                last,
                self.params.k,
            );
        }
        text
    }

    fn write_report(&self, retained_kmers: usize, masked_kmers: usize) -> String {
        format!(
            "# kmerators-wasm report\n\n\
             - k-mer length: {}\n\
             - minimizer length: {}\n\
             - query sequences: {}\n\
             - query k-mers: {}\n\
             - retained k-mers: {}\n\
             - masked k-mers: {}\n\
             - transcriptome filter: {}\n\
             - genome filter: {}\n",
            self.params.k,
            self.params.minimizer_length,
            self.query_items.len(),
            self.occurrences.len(),
            retained_kmers,
            masked_kmers,
            threshold_label(
                self.transcriptome_enabled,
                self.params.max_transcriptome_count
            ),
            threshold_label(self.genome_enabled, self.params.max_genome_count)
        )
    }
}

impl QueryIndex {
    fn new(k: usize) -> Self {
        Self {
            map: HashMap::default(),
            filter: PrefixFilter::new(k),
        }
    }

    fn with_capacity(k: usize, capacity: usize) -> Self {
        Self {
            map: HashMap::with_capacity_and_hasher(capacity, Default::default()),
            filter: PrefixFilter::new(k),
        }
    }

    fn insert(&mut self, minimizer: u64, key: u64, id: usize) {
        self.filter.insert(key);
        self.map
            .entry(minimizer)
            .or_default()
            .entry(key)
            .or_default()
            .push(id);
    }

    fn count_key(&self, minimizer: u64, key: u64, counts: &mut [u32]) {
        if !self.filter.maybe_contains(key) {
            return;
        }
        let Some(partition) = self.map.get(&minimizer) else {
            return;
        };
        if let Some(ids) = partition.get(&key) {
            for &id in ids {
                counts[id] = counts[id].saturating_add(1);
            }
        }
    }
}

impl PrefixFilter {
    fn new(k: usize) -> Self {
        let bits = key_bits(k).min(PREFIX_FILTER_MAX_BITS);
        let word_count = (1usize << bits).div_ceil(64);
        Self {
            k,
            bits,
            words: vec![0; word_count],
        }
    }

    fn insert(&mut self, key: u64) {
        let idx = self.prefix(key);
        self.words[idx / 64] |= 1u64 << (idx % 64);
    }

    fn maybe_contains(&self, key: u64) -> bool {
        let idx = self.prefix(key);
        (self.words[idx / 64] & (1u64 << (idx % 64))) != 0
    }

    fn prefix(&self, key: u64) -> usize {
        let key_bits = key_bits(self.k);
        let value = if key_bits > self.bits {
            key >> (key_bits - self.bits)
        } else {
            key
        };
        (value as usize) & ((1usize << self.bits) - 1)
    }
}

impl FastxStreamScanner {
    fn new(k: usize, minimizer_length: usize) -> Self {
        Self {
            mode: StreamMode::Unknown,
            line: Vec::new(),
            fastq_state: FastqState::Header,
            fastq_seq_len: 0,
            fastq_quality_len: 0,
            roller: KmerRoller::new(k, minimizer_length),
        }
    }

    fn scan(&mut self, bytes: &[u8], key_kind: CountKey, index: &QueryIndex, counts: &mut [u32]) {
        for &byte in bytes {
            if byte == b'\n' {
                self.process_pending_line(key_kind, index, counts);
            } else {
                self.line.push(byte);
            }
        }
    }

    fn finish_source(&mut self, key_kind: CountKey, index: &QueryIndex, counts: &mut [u32]) {
        if !self.line.is_empty() {
            let line = strip_cr(&self.line).to_vec();
            self.line.clear();
            self.process_line(&line, key_kind, index, counts);
        }
        self.mode = StreamMode::Unknown;
        self.fastq_state = FastqState::Header;
        self.fastq_seq_len = 0;
        self.fastq_quality_len = 0;
        self.roller.reset();
    }

    fn process_pending_line(&mut self, key_kind: CountKey, index: &QueryIndex, counts: &mut [u32]) {
        let line = strip_cr(&self.line).to_vec();
        self.line.clear();
        self.process_line(&line, key_kind, index, counts);
    }

    fn process_line(
        &mut self,
        line: &[u8],
        key_kind: CountKey,
        index: &QueryIndex,
        counts: &mut [u32],
    ) {
        let trimmed = trim_ascii(line);
        if trimmed.is_empty() {
            return;
        }

        if self.mode == StreamMode::Unknown {
            self.mode = if trimmed.starts_with(b">") {
                StreamMode::Fasta
            } else if trimmed.starts_with(b"@") {
                StreamMode::Fastq
            } else {
                StreamMode::Fasta
            };
        }

        match self.mode {
            StreamMode::Unknown => {}
            StreamMode::Fasta => self.process_fasta_line(trimmed, key_kind, index, counts),
            StreamMode::Fastq => self.process_fastq_line(line, trimmed, key_kind, index, counts),
        }
    }

    fn process_fasta_line(
        &mut self,
        line: &[u8],
        key_kind: CountKey,
        index: &QueryIndex,
        counts: &mut [u32],
    ) {
        if line.starts_with(b">") {
            self.roller.reset();
            return;
        }
        self.scan_sequence(line, key_kind, index, counts);
    }

    fn process_fastq_line(
        &mut self,
        line: &[u8],
        trimmed: &[u8],
        key_kind: CountKey,
        index: &QueryIndex,
        counts: &mut [u32],
    ) {
        match self.fastq_state {
            FastqState::Header => {
                if trimmed.starts_with(b"@") {
                    self.roller.reset();
                    self.fastq_seq_len = 0;
                    self.fastq_quality_len = 0;
                    self.fastq_state = FastqState::Sequence;
                }
            }
            FastqState::Sequence => {
                if trimmed.starts_with(b"+") {
                    self.fastq_state = FastqState::Quality;
                    self.fastq_quality_len = 0;
                    return;
                }
                self.fastq_seq_len += count_non_whitespace(trimmed);
                self.scan_sequence(trimmed, key_kind, index, counts);
            }
            FastqState::Quality => {
                self.fastq_quality_len += strip_cr(line).len();
                if self.fastq_quality_len >= self.fastq_seq_len {
                    self.fastq_state = FastqState::Header;
                }
            }
        }
    }

    fn scan_sequence(
        &mut self,
        line: &[u8],
        key_kind: CountKey,
        index: &QueryIndex,
        counts: &mut [u32],
    ) {
        for &byte in line {
            if byte.is_ascii_whitespace() {
                continue;
            }
            if let Some(hit) = self.roller.push(byte) {
                let key = match key_kind {
                    CountKey::Forward => hit.forward_key,
                    CountKey::Canonical => hit.canonical_key,
                };
                index.count_key(hit.minimizer, key, counts);
            }
        }
    }
}

impl KmerRoller {
    fn new(k: usize, minimizer_length: usize) -> Self {
        let mode = if k <= MAX_EXACT_K {
            KeyMode::Exact
        } else {
            KeyMode::Hashed
        };
        Self {
            k,
            m: minimizer_length,
            mode,
            valid_len: 0,
            forward: 0,
            reverse: 0,
            exact_mask: exact_mask(k),
            rc_shift: 2 * (k.saturating_sub(1)),
            minimizer_forward: 0,
            minimizer_mask: exact_mask(minimizer_length),
            minimizers: VecDeque::new(),
            window: VecDeque::with_capacity(k),
        }
    }

    fn push(&mut self, byte: u8) -> Option<KmerHit> {
        let bits = match base_bits(byte) {
            Some(bits) => bits,
            None => {
                self.reset();
                return None;
            }
        };
        let base = normalize_base(byte);
        self.valid_len += 1;

        if matches!(self.mode, KeyMode::Exact) {
            self.forward = ((self.forward << 2) | bits) & self.exact_mask;
            self.reverse = (self.reverse >> 2) | ((bits ^ 0b10) << self.rc_shift);
        } else {
            self.window.push_back(base);
            if self.window.len() > self.k {
                self.window.pop_front();
            }
        }

        self.minimizer_forward = ((self.minimizer_forward << 2) | bits) & self.minimizer_mask;
        if self.valid_len >= self.m {
            let mmer_start = self.valid_len - self.m;
            while self
                .minimizers
                .back()
                .is_some_and(|&(_, value)| value > self.minimizer_forward)
            {
                self.minimizers.pop_back();
            }
            self.minimizers
                .push_back((mmer_start, self.minimizer_forward));
        }

        if self.valid_len < self.k {
            return None;
        }

        let kmer_start = self.valid_len - self.k;
        while self
            .minimizers
            .front()
            .is_some_and(|&(pos, _)| pos < kmer_start)
        {
            self.minimizers.pop_front();
        }
        let minimizer = self.minimizers.front().map(|&(_, value)| value)?;
        let (forward_key, canonical_key) = match self.mode {
            KeyMode::Exact => (self.forward, self.forward.min(self.reverse)),
            KeyMode::Hashed => {
                let forward = hash_bases(self.window.iter().copied());
                let reverse = hash_revcomp_bases(self.window.iter().copied().rev());
                (forward, forward.min(reverse))
            }
        };
        Some(KmerHit {
            forward_key,
            canonical_key,
            minimizer,
        })
    }

    fn reset(&mut self) {
        self.valid_len = 0;
        self.forward = 0;
        self.reverse = 0;
        self.minimizer_forward = 0;
        self.minimizers.clear();
        self.window.clear();
    }
}

fn scan_chunk(
    bytes: &[u8],
    key_kind: CountKey,
    index: &QueryIndex,
    counts: &mut [u32],
    scanner: &mut FastxStreamScanner,
) {
    scanner.scan(bytes, key_kind, index, counts);
}

fn parse_query_fastx(bytes: &[u8]) -> Result<Vec<QueryItem>, String> {
    match bytes
        .iter()
        .copied()
        .find(|byte| !byte.is_ascii_whitespace())
    {
        Some(b'>') => parse_query_fasta(bytes),
        Some(b'@') => parse_query_fastq(bytes),
        Some(_) => parse_query_fasta(bytes),
        None => Ok(Vec::new()),
    }
}

fn parse_query_fasta(bytes: &[u8]) -> Result<Vec<QueryItem>, String> {
    let mut items = Vec::new();
    let mut current_id = None::<String>;
    let mut current_seq = Vec::new();
    let mut saw_header = false;

    for raw_line in bytes.split(|&byte| byte == b'\n') {
        let line = raw_line.strip_suffix(b"\r").unwrap_or(raw_line);
        let trimmed = trim_ascii(line);
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix(b">") {
            saw_header = true;
            flush_query_item(&mut items, &mut current_id, &mut current_seq)?;
            current_id = Some(first_token(header).unwrap_or_else(|| {
                let next_id = items.len() + 1;
                format!("query_{next_id}")
            }));
        } else {
            if current_id.is_none() && !saw_header {
                current_id = Some("query_1".to_string());
            }
            for &byte in trimmed {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                current_seq.push(normalize_base(byte));
            }
        }
    }

    flush_query_item(&mut items, &mut current_id, &mut current_seq)?;
    Ok(items)
}

fn parse_query_fastq(bytes: &[u8]) -> Result<Vec<QueryItem>, String> {
    let mut items = Vec::new();
    let mut lines = bytes.split(|&byte| byte == b'\n').peekable();

    loop {
        let Some(header) = next_nonempty_line(&mut lines) else {
            break;
        };
        let header = trim_ascii(strip_cr(header));
        let Some(header) = header.strip_prefix(b"@") else {
            return Err("FASTQ record does not start with @".to_string());
        };
        let id = first_token(header).unwrap_or_else(|| {
            let next_id = items.len() + 1;
            format!("query_{next_id}")
        });

        let mut seq = Vec::new();
        let mut saw_plus = false;
        for raw_line in lines.by_ref() {
            let line = trim_ascii(strip_cr(raw_line));
            if line.starts_with(b"+") {
                saw_plus = true;
                break;
            }
            for &byte in line {
                if !byte.is_ascii_whitespace() {
                    seq.push(normalize_base(byte));
                }
            }
        }
        if !saw_plus {
            return Err(format!("FASTQ record {id} is missing a + separator"));
        }
        if seq.is_empty() {
            return Err(format!("query sequence {id} has no bases"));
        }

        let mut quality_len = 0usize;
        for raw_line in lines.by_ref() {
            quality_len += strip_cr(raw_line).len();
            if quality_len >= seq.len() {
                break;
            }
        }
        if quality_len < seq.len() {
            return Err(format!("FASTQ record {id} has too few quality characters"));
        }

        items.push(QueryItem { id, seq });
    }

    Ok(items)
}

fn flush_query_item(
    items: &mut Vec<QueryItem>,
    current_id: &mut Option<String>,
    current_seq: &mut Vec<u8>,
) -> Result<(), String> {
    let Some(id) = current_id.take() else {
        current_seq.clear();
        return Ok(());
    };
    if current_seq.is_empty() {
        return Err(format!("query sequence {id} has no bases"));
    }
    items.push(QueryItem {
        id,
        seq: std::mem::take(current_seq),
    });
    Ok(())
}

fn build_occurrences(k: usize, minimizer_length: usize, items: &[QueryItem]) -> Vec<Occurrence> {
    let mut occurrences = Vec::new();
    for (item_idx, item) in items.iter().enumerate() {
        let mut roller = KmerRoller::new(k, minimizer_length);
        for (idx, &byte) in item.seq.iter().enumerate() {
            if let Some(hit) = roller.push(byte) {
                let pos = idx + 1 - k;
                occurrences.push(Occurrence {
                    item_idx,
                    pos,
                    forward_key: hit.forward_key,
                    canonical_key: hit.canonical_key,
                    minimizer: hit.minimizer,
                    revcomp_minimizer: revcomp_minimizer(&item.seq[pos..pos + k], minimizer_length)
                        .unwrap_or(hit.minimizer),
                });
            }
        }
    }
    occurrences
}

fn base_bits(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'T' | b't' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

fn normalize_base(base: u8) -> u8 {
    match base {
        b'a'..=b'z' => base - 32,
        _ => base,
    }
}

fn exact_mask(k: usize) -> u64 {
    if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

fn key_bits(k: usize) -> usize {
    if k <= MAX_EXACT_K { 2 * k } else { 64 }
}

fn hash_bases<I>(bases: I) -> u64
where
    I: IntoIterator<Item = u8>,
{
    let mut hash = HASH_OFFSET;
    for base in bases {
        hash ^= u64::from(base_bits(base).unwrap_or(4));
        hash = hash.wrapping_mul(HASH_PRIME);
    }
    hash
}

fn hash_revcomp_bases<I>(bases: I) -> u64
where
    I: IntoIterator<Item = u8>,
{
    let mut hash = HASH_OFFSET;
    for base in bases {
        hash ^= u64::from(base_bits(base).map(|bits| bits ^ 0b10).unwrap_or(4));
        hash = hash.wrapping_mul(HASH_PRIME);
    }
    hash
}

fn revcomp_minimizer(seq: &[u8], m: usize) -> Option<u64> {
    if seq.len() < m {
        return None;
    }
    let mask = exact_mask(m);
    let mut value = 0u64;
    let mut valid_len = 0usize;
    let mut best = None::<u64>;
    for &base in seq.iter().rev() {
        let bits = base_bits(base)? ^ 0b10;
        value = ((value << 2) | bits) & mask;
        valid_len += 1;
        if valid_len >= m && best.is_none_or(|current| value < current) {
            best = Some(value);
        }
    }
    best
}

fn strip_cr(bytes: &[u8]) -> &[u8] {
    bytes.strip_suffix(b"\r").unwrap_or(bytes)
}

fn trim_ascii(bytes: &[u8]) -> &[u8] {
    let start = bytes
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(bytes.len());
    let end = bytes
        .iter()
        .rposition(|byte| !byte.is_ascii_whitespace())
        .map(|idx| idx + 1)
        .unwrap_or(start);
    &bytes[start..end]
}

fn count_non_whitespace(bytes: &[u8]) -> usize {
    bytes
        .iter()
        .filter(|byte| !byte.is_ascii_whitespace())
        .count()
}

fn next_nonempty_line<'a, I>(lines: &mut I) -> Option<&'a [u8]>
where
    I: Iterator<Item = &'a [u8]>,
{
    lines.find(|line| !trim_ascii(strip_cr(line)).is_empty())
}

fn first_token(header: &[u8]) -> Option<String> {
    let token = header
        .split(|byte| byte.is_ascii_whitespace())
        .find(|part| !part.is_empty())?;
    Some(String::from_utf8_lossy(token).into_owned())
}

fn push_fasta_record(text: &mut String, header: &str, seq: &[u8]) {
    text.push('>');
    text.push_str(header);
    text.push('\n');
    text.push_str(&String::from_utf8_lossy(seq));
    text.push('\n');
}

fn count_label(enabled: bool, count: u32) -> String {
    if enabled {
        count.to_string()
    } else {
        "NA".to_string()
    }
}

fn threshold_label(enabled: bool, threshold: u32) -> String {
    if enabled {
        format!("enabled, max count {threshold}")
    } else {
        "skipped".to_string()
    }
}

fn push_contig(
    text: &mut String,
    item: &QueryItem,
    contig_idx: usize,
    start: usize,
    last: usize,
    k: usize,
) {
    let end = last + k;
    push_fasta_record(
        text,
        &format!("{}_contig_{}:{}-{}", item.id, contig_idx, start + 1, end),
        &item.seq[start..end],
    );
}

fn set_last_error(message: impl Into<String>) {
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = message.into().into_bytes();
    });
}

fn clear_last_error() {
    LAST_ERROR.with(|slot| slot.borrow_mut().clear());
}

fn with_session_mut<F>(session_id: u32, action: F) -> i32
where
    F: FnOnce(&mut KmeratorRun) -> Result<(), String>,
{
    clear_last_error();
    let result = SESSIONS.with(|sessions| {
        let mut sessions = sessions.borrow_mut();
        let Some(Some(session)) = session_id
            .checked_sub(1)
            .and_then(|idx| sessions.get_mut(idx as usize))
        else {
            return Err("invalid session id".to_string());
        };
        action(session)
    });

    match result {
        Ok(()) => 0,
        Err(message) => {
            set_last_error(message);
            -1
        }
    }
}

unsafe fn bytes_from_raw<'a>(ptr: *const u8, len: usize) -> Result<&'a [u8], String> {
    if ptr.is_null() && len != 0 {
        return Err("received a null pointer with non-zero length".to_string());
    }
    Ok(unsafe { std::slice::from_raw_parts(ptr, len) })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_alloc(len: usize) -> *mut u8 {
    let mut bytes = Vec::<u8>::with_capacity(len);
    let ptr = bytes.as_mut_ptr();
    std::mem::forget(bytes);
    ptr
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn kmerators_dealloc(ptr: *mut u8, len: usize) {
    if !ptr.is_null() && len != 0 {
        unsafe {
            drop(Vec::from_raw_parts(ptr, 0, len));
        }
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_session_new(
    k: u32,
    max_transcriptome_count: u32,
    max_genome_count: u32,
) -> u32 {
    kmerators_session_new_with_minimizer(k, 0, max_transcriptome_count, max_genome_count)
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_session_new_with_minimizer(
    k: u32,
    minimizer_length: u32,
    max_transcriptome_count: u32,
    max_genome_count: u32,
) -> u32 {
    clear_last_error();
    let run = match KmeratorRun::new(RunParams {
        k: k as usize,
        minimizer_length: minimizer_length as usize,
        max_transcriptome_count,
        max_genome_count,
    }) {
        Ok(run) => run,
        Err(message) => {
            set_last_error(message);
            return 0;
        }
    };

    SESSIONS.with(|sessions| {
        let mut sessions = sessions.borrow_mut();
        if let Some((idx, slot)) = sessions
            .iter_mut()
            .enumerate()
            .find(|(_, slot)| slot.is_none())
        {
            *slot = Some(run);
            (idx + 1) as u32
        } else {
            sessions.push(Some(run));
            sessions.len() as u32
        }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_session_free(session_id: u32) {
    SESSIONS.with(|sessions| {
        if let Some(slot) = session_id.checked_sub(1).and_then(|idx| {
            sessions
                .borrow_mut()
                .get_mut(idx as usize)
                .map(|slot| slot.take())
        }) {
            drop(slot);
        }
    });
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn kmerators_set_query_fasta(
    session_id: u32,
    ptr: *const u8,
    len: usize,
) -> i32 {
    let bytes = match unsafe { bytes_from_raw(ptr, len) } {
        Ok(bytes) => bytes,
        Err(message) => {
            set_last_error(message);
            return -1;
        }
    };
    with_session_mut(session_id, |session| session.set_query_fasta(bytes))
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn kmerators_set_query_fastx(
    session_id: u32,
    ptr: *const u8,
    len: usize,
) -> i32 {
    unsafe { kmerators_set_query_fasta(session_id, ptr, len) }
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_set_transcriptome_enabled(session_id: u32, enabled: u32) -> i32 {
    with_session_mut(session_id, |session| {
        session.set_transcriptome_enabled(enabled != 0);
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_set_genome_enabled(session_id: u32, enabled: u32) -> i32 {
    with_session_mut(session_id, |session| {
        session.set_genome_enabled(enabled != 0);
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn kmerators_scan_transcriptome_chunk(
    session_id: u32,
    ptr: *const u8,
    len: usize,
) -> i32 {
    let bytes = match unsafe { bytes_from_raw(ptr, len) } {
        Ok(bytes) => bytes,
        Err(message) => {
            set_last_error(message);
            return -1;
        }
    };
    with_session_mut(session_id, |session| {
        session.scan_transcriptome_chunk(bytes)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_finish_transcriptome_source(session_id: u32) -> i32 {
    with_session_mut(session_id, |session| {
        session.finish_transcriptome_source();
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn kmerators_scan_genome_chunk(
    session_id: u32,
    ptr: *const u8,
    len: usize,
) -> i32 {
    let bytes = match unsafe { bytes_from_raw(ptr, len) } {
        Ok(bytes) => bytes,
        Err(message) => {
            set_last_error(message);
            return -1;
        }
    };
    with_session_mut(session_id, |session| session.scan_genome_chunk(bytes))
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_finish_genome_source(session_id: u32) -> i32 {
    with_session_mut(session_id, |session| {
        session.finish_genome_source();
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_finish(session_id: u32) -> i32 {
    with_session_mut(session_id, |session| session.finish().map(|_| ()))
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_result_ptr(session_id: u32) -> *const u8 {
    SESSIONS.with(|sessions| {
        let sessions = sessions.borrow();
        session_id
            .checked_sub(1)
            .and_then(|idx| sessions.get(idx as usize))
            .and_then(|slot| slot.as_ref())
            .map(|session| session.last_json().as_ptr())
            .unwrap_or(std::ptr::null())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_result_len(session_id: u32) -> usize {
    SESSIONS.with(|sessions| {
        let sessions = sessions.borrow();
        session_id
            .checked_sub(1)
            .and_then(|idx| sessions.get(idx as usize))
            .and_then(|slot| slot.as_ref())
            .map(|session| session.last_json().len())
            .unwrap_or(0)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_last_error_ptr() -> *const u8 {
    LAST_ERROR.with(|slot| slot.borrow().as_ptr())
}

#[unsafe(no_mangle)]
pub extern "C" fn kmerators_last_error_len() -> usize {
    LAST_ERROR.with(|slot| slot.borrow().len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_tiny_with_chunks(chunks: &[&[u8]]) -> RunOutput {
        let mut run = KmeratorRun::new(RunParams {
            k: 5,
            minimizer_length: 0,
            max_transcriptome_count: 0,
            max_genome_count: 1,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nACGTACCCC\n").unwrap();
        run.scan_transcriptome_chunk(b">tr1\nACGTA\n").unwrap();
        for chunk in chunks {
            run.scan_genome_chunk(chunk).unwrap();
        }
        run.finish().unwrap();
        serde_json::from_slice(run.last_json()).unwrap()
    }

    #[test]
    fn streaming_scan_matches_unsplit_reference() {
        let full = run_tiny_with_chunks(&[b">chr1\nACGTACCCC\n"]);
        let split = run_tiny_with_chunks(&[b">chr1\nAC", b"GTAC", b"CCC\n"]);

        assert_eq!(full.retained_kmers, split.retained_kmers);
        assert_eq!(full.masked_kmers, split.masked_kmers);
        assert_eq!(full.files.contigs_fa, split.files.contigs_fa);
    }

    #[test]
    fn invalid_bases_break_kmer_windows() {
        let mut run = KmeratorRun::new(RunParams {
            k: 3,
            minimizer_length: 0,
            max_transcriptome_count: 99,
            max_genome_count: 99,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nACG\n").unwrap();
        run.scan_genome_chunk(b">chr1\nACN").unwrap();
        run.scan_genome_chunk(b"G\n").unwrap();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.genome_hits, Some(0));
    }

    #[test]
    fn headers_break_kmer_windows() {
        let mut run = KmeratorRun::new(RunParams {
            k: 4,
            minimizer_length: 0,
            max_transcriptome_count: 99,
            max_genome_count: 99,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nACGT\n").unwrap();
        run.scan_genome_chunk(b">chr1\nAC\n>chr2\nGT\n").unwrap();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.genome_hits, Some(0));
    }

    #[test]
    fn reverse_complement_counts_for_genome() {
        let mut run = KmeratorRun::new(RunParams {
            k: 5,
            minimizer_length: 0,
            max_transcriptome_count: 99,
            max_genome_count: 99,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nAACGT\n").unwrap();
        run.scan_genome_chunk(b">chr1\nACGTT\n").unwrap();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.genome_hits, Some(1));
    }

    #[test]
    fn retained_kmers_merge_into_contigs() {
        let mut run = KmeratorRun::new(RunParams {
            k: 3,
            minimizer_length: 0,
            max_transcriptome_count: 0,
            max_genome_count: 10,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nACGTA\n").unwrap();
        run.scan_genome_chunk(b">chr1\nACGTA\n").unwrap();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert!(
            output
                .files
                .contigs_fa
                .contains(">q1_contig_1:1-5\nACGTA\n")
        );
    }

    #[test]
    fn skipped_references_do_not_filter_or_fake_counts() {
        let mut run = KmeratorRun::new(RunParams {
            k: 3,
            minimizer_length: 0,
            max_transcriptome_count: 0,
            max_genome_count: 0,
        })
        .unwrap();
        run.set_query_fasta(b">q1\nACGTA\n").unwrap();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.retained_kmers, 3);
        assert_eq!(output.transcriptome_hits, None);
        assert_eq!(output.genome_hits, None);
        assert!(
            output
                .files
                .report_md
                .contains("transcriptome filter: skipped")
        );
        assert!(output.files.report_md.contains("genome filter: skipped"));
    }

    #[test]
    fn fastq_query_and_reference_ignore_quality_lines() {
        let mut run = KmeratorRun::new(RunParams {
            k: 5,
            minimizer_length: 0,
            max_transcriptome_count: 99,
            max_genome_count: 99,
        })
        .unwrap();
        run.set_query_fasta(b"@q1\nACGTA\n+\n!!!!!\n").unwrap();
        run.scan_genome_chunk(b"@chr1\nACGTA\n+\nAAAAA\n").unwrap();
        run.finish_genome_source();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.query_kmers, 1);
        assert_eq!(output.genome_hits, Some(1));
    }

    #[test]
    fn long_kmers_use_hashed_keys() {
        let query = b"ACGTTGCAACGTGGTACCTTAGGCTAACCGTATGCCGTAACCTGG";
        let mut query_fasta = b">q1\n".to_vec();
        query_fasta.extend_from_slice(query);
        query_fasta.push(b'\n');

        let mut genome_fasta = b">chr1\n".to_vec();
        genome_fasta.extend_from_slice(query);

        let mut run = KmeratorRun::new(RunParams {
            k: 41,
            minimizer_length: 9,
            max_transcriptome_count: 99,
            max_genome_count: 99,
        })
        .unwrap();
        run.set_query_fasta(&query_fasta).unwrap();
        run.scan_genome_chunk(&genome_fasta).unwrap();
        run.finish_genome_source();
        run.finish().unwrap();
        let output: RunOutput = serde_json::from_slice(run.last_json()).unwrap();

        assert_eq!(output.query_kmers, query.len() - 41 + 1);
        assert_eq!(output.genome_hits, Some((query.len() - 41 + 1) as u64));
        assert!(output.files.report_md.contains("minimizer length: 9"));
    }
}
