use crackle_kit::tracing::{Level, event};

// 1. Make the Enum PRIVATE (no 'pub')
// The user cannot see this or type 'Strategy::Simd'
#[derive(Clone, Copy)]
enum Strategy {
    Normal,
    Simd,
}

// 2. Make a Public Struct that WRAPS the Enum
// This is the only thing the user sees.
#[derive(Clone)]
pub struct GCATCounter {
    inner: Strategy,
}
impl GCATCounter {
    // 3. The Constructor handles the safety check ONCE.
    pub fn new() -> Self {
        let strategy = if is_x86_feature_detected!("avx2") {
            event!(Level::DEBUG, "avx2 enabled. use simd mode.");
            Strategy::Simd
            // Strategy::Normal
        } else {
            event!(Level::DEBUG, "avx2 not found. use normal mode.");
            Strategy::Normal
        };

        Self { inner: strategy }
    }

    // 4. The public method just delegates
    pub fn count(&self, ref_seq: &[u8]) -> (i32, i32) {
        match self.inner {
            Strategy::Normal => count_gc_and_at(ref_seq),
            Strategy::Simd => {
                // SAFETY: We checked AVX2 inside 'new()', so this is safe.
                #[cfg(all(target_arch = "x86_64"))]
                {
                    unsafe { count_gc_simd_avx2(ref_seq) }
                }

                #[cfg(not(all(target_arch = "x86_64")))]
                {
                    count_gc_and_at(ref_seq)
                }
            }
        }
    }
}

// 5. Default trait just calls new()
impl Default for GCATCounter {
    fn default() -> Self {
        Self::new()
    }
}

// pub fn calculate_gc_at(bytes: &[u8], lookup_table: &[u8; 256]) -> (i32, i32) {
//     #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
//     {
//         // If compiled with AVX2 support, use the SIMD monster
//         count_gc_simd_avx2(bytes, lookup_table)
//     }

//     #[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
//     {
//         // Otherwise, use the scalar fallback
//         count_gc_and_at(bytes, lookup_table)
//     }
// }

#[target_feature(enable = "avx2")]
#[cfg(all(target_arch = "x86_64"))]
pub(crate) fn count_gc_simd_avx2(bytes: &[u8]) -> (i32, i32) {
    use std::arch::x86_64::*;
    let len = bytes.len();
    let num_iter = len / 32;
    let chunks_len = num_iter * 32;
    let mut gc_count = 0;
    let mut at_count = 0;

    unsafe {
        let lowercase_bit = _mm256_set1_epi8(0x20); // 0b0010_0000
        let needle_g = _mm256_set1_epi8(b'g' as i8);
        let needle_c = _mm256_set1_epi8(b'c' as i8);
        let needle_a = _mm256_set1_epi8(b'a' as i8);
        let needle_t = _mm256_set1_epi8(b't' as i8);

        for i in 0..num_iter {
            let chunk = _mm256_loadu_si256(bytes.as_ptr().add(i * 32) as *const __m256i);

            // --- THE FAIR NORMALIZATION ---
            // Forced-set the 6th bit.
            // Now 'G' and 'g' both become 'g'. 'A' and 'a' both become 'a'.
            let normalized = _mm256_or_si256(chunk, lowercase_bit);

            // One check covers BOTH cases
            let gc_mask = _mm256_or_si256(
                _mm256_cmpeq_epi8(normalized, needle_g),
                _mm256_cmpeq_epi8(normalized, needle_c),
            );
            let at_mask = _mm256_or_si256(
                _mm256_cmpeq_epi8(normalized, needle_a),
                _mm256_cmpeq_epi8(normalized, needle_t),
            );

            gc_count += _mm256_movemask_epi8(gc_mask).count_ones() as usize;
            at_count += _mm256_movemask_epi8(at_mask).count_ones() as usize;
        }
    }

    let (t_gc, t_at) = count_gc_and_at(&bytes[chunks_len..]);
    (gc_count as i32 + t_gc, at_count as i32 + t_at)
}

const GC_TABLE: [u8; 256] = {
        let mut table = [0; 256];
        // Mark GC as 1, AT as 2, others 0
        table[b'G' as usize] = 1;
        table[b'g' as usize] = 1;
        table[b'C' as usize] = 1;
        table[b'c' as usize] = 1;
        table[b'A' as usize] = 2;
        table[b'a' as usize] = 2;
        table[b'T' as usize] = 2;
        table[b't' as usize] = 2;
        table
    };

pub(crate) fn count_gc_and_at(ref_seq: &[u8]) -> (i32, i32) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for &b in ref_seq {
        match GC_TABLE[b as usize] {
            1 => gc_cnt += 1,
            2 => at_cnt += 1,
            _ => {}
        }
    }

    (gc_cnt, at_cnt)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to generate random DNA-like bytes
    fn generate_random_dna(len: usize) -> Vec<u8> {
        use rand::Rng; // Add 'rand' to dev-dependencies in Cargo.toml
        let mut rng = rand::thread_rng();
        let bases = b"ACGTacgtN"; // Include mixed case and junk
        (0..len)
            .map(|_| bases[rng.gen_range(0..bases.len())])
            .collect()
    }

    
    #[test]
    fn test_equivalence_short() {
        // Case 1: Short string (Scalar fallback only)
        let data = b"ACGTacgt"; // 8 bytes
        
        // Only run AVX2 test if CPU supports it
        if is_x86_feature_detected!("avx2") {
            let simd_res = unsafe { count_gc_simd_avx2(data) };
            let scalar_res = count_gc_and_at(data);
            assert_eq!(simd_res, scalar_res, "Failed on short string");
        }
    }

    #[test]
    fn test_equivalence_exact_chunk() {
        // Case 2: Exactly 32 bytes (1 SIMD loop, 0 scalar tail)
        let data = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; 
        
        if is_x86_feature_detected!("avx2") {
            let simd_res = unsafe { count_gc_simd_avx2(data) };
            let scalar_res = count_gc_and_at(data);
            assert_eq!(simd_res, scalar_res, "Failed on exact 32-byte chunk");
        }
    }

    #[test]
    fn test_equivalence_chunk_plus_one() {
        // Case 3: 33 bytes (1 SIMD loop, 1 scalar tail)
        let data = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG"; 
        
        if is_x86_feature_detected!("avx2") {
            let simd_res = unsafe { count_gc_simd_avx2(data) };
            let scalar_res = count_gc_and_at(data);
            assert_eq!(simd_res, scalar_res, "Failed on 33 bytes");
        }
    }

    #[test]
    fn test_equivalence_fuzzing() {
        if !is_x86_feature_detected!("avx2") {
            println!("Skipping AVX2 fuzzing test on non-AVX2 machine");
            return;
        }

        // Run 100 random tests with varying lengths
        for _ in 0..100 {
            use rand::Rng;
            let mut rng = rand::thread_rng();
            // Random length between 0 and 2000 to test many chunk boundaries
            let len = rng.gen_range(0..2000); 
            let data = generate_random_dna(len);

            let simd_res = unsafe { count_gc_simd_avx2(&data) };
            let scalar_res = count_gc_and_at(&data);

            assert_eq!(
                simd_res, 
                scalar_res, 
                "Mismatch found! Length: {}, Data: {:?}", 
                len, 
                String::from_utf8_lossy(&data)
            );
        }
    }
}