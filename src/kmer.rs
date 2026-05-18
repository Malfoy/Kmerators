use ahash::AHasher;
use std::hash::Hasher;

pub const MAX_EXACT_U64_K: usize = 31;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KeyMode {
    ExactU64,
    HashedU64,
}

impl KeyMode {
    pub fn for_k(k: usize) -> Self {
        if k <= MAX_EXACT_U64_K {
            Self::ExactU64
        } else {
            Self::HashedU64
        }
    }
}

pub fn base_bits(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'T' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

pub fn encode_forward_exact(seq: &[u8]) -> Option<u64> {
    let mut value = 0u64;
    for &base in seq {
        value = (value << 2) | base_bits(base)?;
    }
    Some(value)
}

pub fn encode_revcomp_exact(seq: &[u8]) -> Option<u64> {
    let mut value = 0u64;
    for &base in seq.iter().rev() {
        value = (value << 2) | (base_bits(base)? ^ 0b10);
    }
    Some(value)
}

pub fn forward_key(seq: &[u8], mode: KeyMode) -> Option<u64> {
    match mode {
        KeyMode::ExactU64 => encode_forward_exact(seq),
        KeyMode::HashedU64 => Some(hash_bytes(seq)),
    }
}

pub fn canonical_key(seq: &[u8], mode: KeyMode) -> Option<u64> {
    match mode {
        KeyMode::ExactU64 => {
            let fwd = encode_forward_exact(seq)?;
            let rc = encode_revcomp_exact(seq)?;
            Some(fwd.min(rc))
        }
        KeyMode::HashedU64 => {
            let fwd = hash_bytes(seq);
            let rc = hash_revcomp(seq)?;
            Some(fwd.min(rc))
        }
    }
}

pub fn revcomp(seq: &[u8]) -> Option<Vec<u8>> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => Some(b'T'),
            b'C' => Some(b'G'),
            b'T' => Some(b'A'),
            b'G' => Some(b'C'),
            _ => None,
        })
        .collect()
}

pub fn hash_bytes(seq: &[u8]) -> u64 {
    let mut hasher = AHasher::default();
    hasher.write_usize(seq.len());
    for &base in seq {
        hasher.write_u8(base);
    }
    hasher.finish()
}

pub fn hash_revcomp(seq: &[u8]) -> Option<u64> {
    let mut hasher = AHasher::default();
    hasher.write_usize(seq.len());
    for &base in seq.iter().rev() {
        let comp = match base {
            b'A' => b'T',
            b'C' => b'G',
            b'T' => b'A',
            b'G' => b'C',
            _ => return None,
        };
        hasher.write_u8(comp);
    }
    Some(hasher.finish())
}

pub fn valid_chunks(seq: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    struct Chunks<'a> {
        seq: &'a [u8],
        pos: usize,
    }

    impl<'a> Iterator for Chunks<'a> {
        type Item = (usize, &'a [u8]);

        fn next(&mut self) -> Option<Self::Item> {
            while self.pos < self.seq.len() && base_bits(self.seq[self.pos]).is_none() {
                self.pos += 1;
            }
            if self.pos >= self.seq.len() {
                return None;
            }
            let start = self.pos;
            while self.pos < self.seq.len() && base_bits(self.seq[self.pos]).is_some() {
                self.pos += 1;
            }
            Some((start, &self.seq[start..self.pos]))
        }
    }

    Chunks { seq, pos: 0 }
}

pub fn minimizer_values(seq: &[u8], m: usize) -> Vec<u64> {
    if seq.len() < m {
        return Vec::new();
    }
    let mask = if m == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * m)) - 1
    };
    let mut value = 0u64;
    let mut values = Vec::with_capacity(seq.len() - m + 1);
    for (idx, &base) in seq.iter().enumerate() {
        value = ((value << 2) | base_bits(base).expect("minimizer input must be ACTG")) & mask;
        if idx + 1 >= m {
            values.push(value);
        }
    }
    values
}

pub fn window_minima(values: &[u64], window: usize) -> Vec<u64> {
    use std::collections::VecDeque;

    if window == 0 || values.len() < window {
        return Vec::new();
    }

    let mut deque: VecDeque<usize> = VecDeque::new();
    let mut minima = Vec::with_capacity(values.len() - window + 1);

    for idx in 0..values.len() {
        while let Some(&back) = deque.back() {
            if values[back] <= values[idx] {
                break;
            }
            deque.pop_back();
        }
        deque.push_back(idx);

        if let Some(&front) = deque.front()
            && front + window <= idx
        {
            deque.pop_front();
        }

        if idx + 1 >= window {
            let front = *deque.front().expect("deque is not empty");
            minima.push(values[front]);
        }
    }

    minima
}

pub fn kmer_keys_for_chunk(seq: &[u8], k: usize, mode: KeyMode) -> Vec<(u64, u64)> {
    if seq.len() < k {
        return Vec::new();
    }

    match mode {
        KeyMode::ExactU64 => exact_keys_for_chunk(seq, k),
        KeyMode::HashedU64 => (0..=seq.len() - k)
            .filter_map(|pos| {
                let kmer = &seq[pos..pos + k];
                Some((forward_key(kmer, mode)?, canonical_key(kmer, mode)?))
            })
            .collect(),
    }
}

fn exact_keys_for_chunk(seq: &[u8], k: usize) -> Vec<(u64, u64)> {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_shift = 2 * (k - 1);
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut keys = Vec::with_capacity(seq.len() - k + 1);

    for (idx, &base) in seq.iter().enumerate() {
        let bits = base_bits(base).expect("exact input must be ACTG");
        fwd = ((fwd << 2) | bits) & mask;
        rev = (rev >> 2) | ((bits ^ 0b10) << rc_shift);
        if idx + 1 >= k {
            keys.push((fwd, fwd.min(rev)));
        }
    }

    keys
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_canonical_matches_revcomp() {
        let seq = b"ACGTT";
        let rc = revcomp(seq).unwrap();
        assert_eq!(
            canonical_key(seq, KeyMode::ExactU64),
            canonical_key(&rc, KeyMode::ExactU64)
        );
    }

    #[test]
    fn finds_valid_chunks() {
        let chunks = valid_chunks(b"ACNNtgAC")
            .map(|(pos, chunk)| (pos, chunk.to_vec()))
            .collect::<Vec<_>>();
        assert_eq!(chunks, vec![(0, b"AC".to_vec()), (6, b"AC".to_vec())]);
    }

    #[test]
    fn sliding_window_minima() {
        assert_eq!(window_minima(&[5, 4, 7, 3, 6], 3), vec![4, 3, 3]);
    }

    #[test]
    fn key_mode_uses_exact_keys_through_k31() {
        assert_eq!(KeyMode::for_k(31), KeyMode::ExactU64);
        assert_eq!(KeyMode::for_k(32), KeyMode::HashedU64);
    }

    #[test]
    fn invalid_bases_reject_exact_keys_and_revcomp() {
        assert_eq!(encode_forward_exact(b"ACN"), None);
        assert_eq!(canonical_key(b"ACN", KeyMode::ExactU64), None);
        assert_eq!(revcomp(b"ACN"), None);
    }

    #[test]
    fn minimizer_values_use_two_bit_rolling_encoding() {
        assert_eq!(minimizer_values(b"ACTG", 2), vec![1, 6, 11]);
    }

    #[test]
    fn kmer_keys_include_forward_and_canonical_values() {
        assert_eq!(
            kmer_keys_for_chunk(b"ACGT", 3, KeyMode::ExactU64),
            vec![(7, 7), (30, 7)]
        );
    }

    #[test]
    fn hashed_canonical_matches_reverse_complement_for_long_kmers() {
        let seq = b"ACTGACTGACTGACTGACTGACTGACTGACTG";
        let rc = revcomp(seq).unwrap();

        assert_eq!(
            canonical_key(seq, KeyMode::HashedU64),
            canonical_key(&rc, KeyMode::HashedU64)
        );
    }
}
