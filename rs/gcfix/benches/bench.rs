use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rand::Rng;

// --- 1. Original (to_ascii_lowercase) ---
fn count_gc_conversion(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for b in seq {
        match b.to_ascii_lowercase() {
            b'g' | b'c' => gc_cnt += 1,
            b'a' | b't' => at_cnt += 1,
            _ => continue,
        }
    }
    (gc_cnt, at_cnt)
}

// --- 2. Optimized (Direct Match) ---
fn count_gc_direct_match(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => gc_cnt += 1,
            b'A' | b'a' | b'T' | b't' => at_cnt += 1,
            _ => continue,
        }
    }
    (gc_cnt, at_cnt)
}

// --- 3. Lookup Table ---
const GC_TABLE: [u8; 256] = {
    let mut table = [0; 256];
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

fn count_gc_lookup(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for &b in seq {
        match GC_TABLE[b as usize] {
            1 => gc_cnt += 1,
            2 => at_cnt += 1,
            _ => {}
        }
    }
    (gc_cnt, at_cnt)
}

// --- 4. Memchr (Crate) ---
fn count_gc_memchr(seq: &[u8]) -> (usize, usize) {
    let gc = memchr::memchr2_iter(b'G', b'g', seq).count()
        + memchr::memchr2_iter(b'C', b'c', seq).count();
    let at = memchr::memchr2_iter(b'A', b'a', seq).count()
        + memchr::memchr2_iter(b'T', b't', seq).count();

    (gc, at)
}

#[cfg(all(target_arch = "x86_64", target_feature="avx"))]
unsafe fn count_gc_simd_avx2(bytes: &[u8]) -> (usize, usize) {
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

    let (t_gc, t_at) = count_gc_direct_match(&bytes[chunks_len..]);
    (gc_count + t_gc, at_count + t_at)
}

fn criterion_benchmark_len_10000(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let bases = b"AGCTagctNn";
    let sequence: Vec<u8> = (0..10_000)
        .map(|_| bases[rng.gen_range(0..bases.len())])
        .collect();

    let mut group = c.benchmark_group("GC Counting Strategies len-10000");

    group.bench_function("1_to_ascii_lowercase", |b| {
        b.iter(|| count_gc_conversion(black_box(&sequence)))
    });
    group.bench_function("2_direct_match_bytes", |b| {
        b.iter(|| count_gc_direct_match(black_box(&sequence)))
    });
    group.bench_function("3_lookup_table", |b| {
        b.iter(|| count_gc_lookup(black_box(&sequence)))
    });
    group.bench_function("4_memchr", |b| {
        b.iter(|| count_gc_memchr(black_box(&sequence)))
    });

    #[cfg(target_arch = "x86_64")]
    group.bench_function("5_avx2_manual", |b| {
        b.iter(|| unsafe { count_gc_simd_avx2(black_box(&sequence)) })
    });

    group.finish();
}

fn criterion_benchmark_len_167(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let bases = b"AGCTagctNn";
    let sequence: Vec<u8> = (0..167)
        .map(|_| bases[rng.gen_range(0..bases.len())])
        .collect();

    let mut group = c.benchmark_group("GC Counting Strategies len-167");

    group.bench_function("1_to_ascii_lowercase", |b| {
        b.iter(|| count_gc_conversion(black_box(&sequence)))
    });
    group.bench_function("2_direct_match_bytes", |b| {
        b.iter(|| count_gc_direct_match(black_box(&sequence)))
    });
    group.bench_function("3_lookup_table", |b| {
        b.iter(|| count_gc_lookup(black_box(&sequence)))
    });
    group.bench_function("4_memchr", |b| {
        b.iter(|| count_gc_memchr(black_box(&sequence)))
    });

    #[cfg(target_arch = "x86_64")]
    group.bench_function("5_avx2_manual", |b| {
        b.iter(|| unsafe { count_gc_simd_avx2(black_box(&sequence)) })
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark_len_10000, criterion_benchmark_len_167);
criterion_main!(benches);
