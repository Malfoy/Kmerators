use serde::{Deserialize, Serialize};
use std::cell::RefCell;

use rustc_hash::FxHashMap as HashMap;

const MAX_EXACT_K: usize = 31;
const PREFIX_FILTER_MAX_BITS: usize = 20;

thread_local! {
    static SESSIONS: RefCell<Vec<Option<KmeratorRun>>> = const { RefCell::new(Vec::new()) };
    static LAST_ERROR: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[derive(Debug, Clone, Copy)]
pub struct RunParams {
    pub k: usize,
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
    transcriptome_scanner: FastaStreamScanner,
    genome_scanner: FastaStreamScanner,
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
}

#[derive(Debug, Clone)]
struct QueryIndex {
    map: HashMap<u64, Vec<usize>>,
    filter: PrefixFilter,
}

#[derive(Debug, Clone)]
struct PrefixFilter {
    k: usize,
    bits: usize,
    words: Vec<u64>,
}

#[derive(Debug, Clone, Copy)]
enum CountKey {
    Forward,
    Canonical,
}

#[derive(Debug)]
struct FastaStreamScanner {
    in_header: bool,
    line_start: bool,
    roller: RollingKmer,
}

#[derive(Debug)]
struct RollingKmer {
    k: usize,
    valid_len: usize,
    forward: u64,
    reverse: u64,
    mask: u64,
    rc_shift: usize,
}

#[derive(Debug, Deserialize, Serialize)]
struct RunOutput {
    kmer_length: usize,
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
        if self.k > MAX_EXACT_K {
            return Err(format!(
                "this browser prototype supports exact k-mers up to {MAX_EXACT_K}"
            ));
        }
        Ok(self)
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
            transcriptome_scanner: FastaStreamScanner::new(params.k),
            genome_scanner: FastaStreamScanner::new(params.k),
            last_json: Vec::new(),
        })
    }

    pub fn set_query_fasta(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.query_items = parse_query_fasta(bytes)?;
        if self.query_items.is_empty() {
            return Err("query FASTA contains no usable sequence".to_string());
        }

        self.occurrences = build_occurrences(self.params.k, &self.query_items);
        if self.occurrences.is_empty() {
            return Err("query FASTA produced no valid k-mers".to_string());
        }

        self.transcriptome_index = QueryIndex::with_capacity(self.params.k, self.occurrences.len());
        self.genome_index = QueryIndex::with_capacity(self.params.k, self.occurrences.len());
        for (idx, occ) in self.occurrences.iter().enumerate() {
            self.transcriptome_index.insert(occ.forward_key, idx);
            self.genome_index.insert(occ.canonical_key, idx);
        }

        self.transcriptome_counts = vec![0; self.occurrences.len()];
        self.genome_counts = vec![0; self.occurrences.len()];
        self.transcriptome_scanner = FastaStreamScanner::new(self.params.k);
        self.genome_scanner = FastaStreamScanner::new(self.params.k);
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
            Err("load a query FASTA before scanning references".to_string())
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
             - query sequences: {}\n\
             - query k-mers: {}\n\
             - retained k-mers: {}\n\
             - masked k-mers: {}\n\
             - transcriptome filter: {}\n\
             - genome filter: {}\n",
            self.params.k,
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

    fn insert(&mut self, key: u64, id: usize) {
        self.filter.insert(key);
        self.map.entry(key).or_default().push(id);
    }

    fn count_key(&self, key: u64, counts: &mut [u32]) {
        if !self.filter.maybe_contains(key) {
            return;
        }
        if let Some(ids) = self.map.get(&key) {
            for &id in ids {
                counts[id] = counts[id].saturating_add(1);
            }
        }
    }
}

impl PrefixFilter {
    fn new(k: usize) -> Self {
        let bits = (2 * k).min(PREFIX_FILTER_MAX_BITS);
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
        let key_bits = 2 * self.k;
        let value = if key_bits > self.bits {
            key >> (key_bits - self.bits)
        } else {
            key
        };
        (value as usize) & ((1usize << self.bits) - 1)
    }
}

impl FastaStreamScanner {
    fn new(k: usize) -> Self {
        Self {
            in_header: false,
            line_start: true,
            roller: RollingKmer::new(k),
        }
    }

    fn scan(&mut self, bytes: &[u8], key_kind: CountKey, index: &QueryIndex, counts: &mut [u32]) {
        for &byte in bytes {
            match byte {
                b'\n' => {
                    if self.in_header {
                        self.in_header = false;
                    }
                    self.line_start = true;
                }
                b'\r' | b' ' | b'\t' => {}
                b'>' if self.line_start => {
                    self.in_header = true;
                    self.line_start = false;
                    self.roller.reset();
                }
                _ if self.in_header => {
                    self.line_start = false;
                }
                _ => {
                    self.line_start = false;
                    if let Some(bits) = base_bits(byte) {
                        if let Some((forward, canonical)) = self.roller.push(bits) {
                            let key = match key_kind {
                                CountKey::Forward => forward,
                                CountKey::Canonical => canonical,
                            };
                            index.count_key(key, counts);
                        }
                    } else {
                        self.roller.reset();
                    }
                }
            }
        }
    }
}

impl RollingKmer {
    fn new(k: usize) -> Self {
        Self {
            k,
            valid_len: 0,
            forward: 0,
            reverse: 0,
            mask: (1u64 << (2 * k)) - 1,
            rc_shift: 2 * (k - 1),
        }
    }

    fn push(&mut self, bits: u64) -> Option<(u64, u64)> {
        self.forward = ((self.forward << 2) | bits) & self.mask;
        self.reverse = (self.reverse >> 2) | ((bits ^ 0b10) << self.rc_shift);
        self.valid_len = (self.valid_len + 1).min(self.k);
        (self.valid_len == self.k).then_some((self.forward, self.forward.min(self.reverse)))
    }

    fn reset(&mut self) {
        self.valid_len = 0;
        self.forward = 0;
        self.reverse = 0;
    }
}

fn scan_chunk(
    bytes: &[u8],
    key_kind: CountKey,
    index: &QueryIndex,
    counts: &mut [u32],
    scanner: &mut FastaStreamScanner,
) {
    scanner.scan(bytes, key_kind, index, counts);
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

fn build_occurrences(k: usize, items: &[QueryItem]) -> Vec<Occurrence> {
    let mut occurrences = Vec::new();
    for (item_idx, item) in items.iter().enumerate() {
        let mut roller = RollingKmer::new(k);
        for (idx, &byte) in item.seq.iter().enumerate() {
            let Some(bits) = base_bits(byte) else {
                roller.reset();
                continue;
            };
            if let Some((forward, canonical)) = roller.push(bits) {
                occurrences.push(Occurrence {
                    item_idx,
                    pos: idx + 1 - k,
                    forward_key: forward,
                    canonical_key: canonical,
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
    clear_last_error();
    let run = match KmeratorRun::new(RunParams {
        k: k as usize,
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
        assert!(output.files.report_md.contains("transcriptome filter: skipped"));
        assert!(output.files.report_md.contains("genome filter: skipped"));
    }
}
