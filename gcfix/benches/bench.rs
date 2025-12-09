use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;

// Method 1: The Original (to_ascii_lowercase)
fn count_gc_conversion(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for b in seq {
        // The conversion happens here every iteration
        match b.to_ascii_lowercase() {
            b'g' | b'c' => {
                gc_cnt += 1;
            }
            b'a' | b't' => {
                at_cnt += 1;
            }
            _ => continue,
        }
    }
    (gc_cnt, at_cnt)
}

// Method 2: The Optimized (Direct Match)
fn count_gc_direct_match(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for &b in seq {
        // No conversion, just a slightly larger jump table
        match b {
            b'G' | b'g' | b'C' | b'c' => {
                gc_cnt += 1;
            }
            b'A' | b'a' | b'T' | b't' => {
                at_cnt += 1;
            }
            _ => continue,
        }
    }
    (gc_cnt, at_cnt)
}

// Method 3: The "Lookup Table" (Simulated Hardcore Optimization)
// This is often used in high-performance tools (like samtools C code)
const GC_TABLE: [u8; 256] = {
    let mut table = [0; 256];
    // Mark GC as 1, AT as 2, others 0
    table[b'G' as usize] = 1; table[b'g' as usize] = 1;
    table[b'C' as usize] = 1; table[b'c' as usize] = 1;
    table[b'A' as usize] = 2; table[b'a' as usize] = 2;
    table[b'T' as usize] = 2; table[b't' as usize] = 2;
    table
};

fn count_gc_lookup(seq: &[u8]) -> (usize, usize) {
    let mut gc_cnt = 0;
    let mut at_cnt = 0;
    for &b in seq {
        // Array indexing is usually faster than branching logic
        match GC_TABLE[b as usize] {
            1 => gc_cnt += 1,
            2 => at_cnt += 1,
            _ => {}
        }
    }
    (gc_cnt, at_cnt)
}

fn criterion_benchmark(c: &mut Criterion) {
    // 1. Generate Random DNA (10kb) with mixed case
    let mut rng = rand::thread_rng();
    let bases = b"AGCTagctNn";
    let sequence: Vec<u8> = (0..10_000)
        .map(|_| bases[rng.gen_range(0..bases.len())])
        .collect();

    let mut group = c.benchmark_group("GC Counting Strategies");

    // Benchmark 1: Conversion
    group.bench_function("to_ascii_lowercase", |b| {
        b.iter(|| count_gc_conversion(black_box(&sequence)))
    });

    // Benchmark 2: Direct Match
    group.bench_function("direct_match_bytes", |b| {
        b.iter(|| count_gc_direct_match(black_box(&sequence)))
    });

    // Benchmark 3: Lookup Table
    group.bench_function("lookup_table", |b| {
        b.iter(|| count_gc_lookup(black_box(&sequence)))
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);