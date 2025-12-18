use std::{
    borrow::Cow,
    cell::RefCell,
    collections::{HashMap, HashSet},
    fs::File,
    i32,
    io::{BufRead, BufReader},
    ops::AddAssign,
    os::unix::thread,
    path::{Path, PathBuf},
};

use std::fmt::Write;

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::locus::GenomeRegion,
    rust_htslib::{
        bam::{IndexedReader, Read, Record, ext::BamRecordExtensions as _, pileup::PileupOption},
        faidx,
    },
    tracing::{Level, event},
};
use ndarray::{Array, Array1, Array2, ArrayBase, Dim, OwnedRepr};
use rayon::{
    ThreadPoolBuilder,
    iter::{
        IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
    },
};

#[derive(Clone, Copy, Debug)]
enum ContigNameFormat {
    WithChr,
    WithoutChr,
}

impl ContigNameFormat {
    fn contig_name_for_ref<'a>(&self, contig: &'a str) -> Cow<'a, str> {
        match self {
            ContigNameFormat::WithChr => {
                if contig.starts_with("chr") {
                    contig.into()
                } else {
                    format!("chr{contig}").into()
                }
            }
            ContigNameFormat::WithoutChr => {
                if contig.starts_with("chr") {
                    contig.get(3..).unwrap().into()
                } else {
                    contig.into()
                }
            }
        }
    }
}

fn make_ref_seq_map(
    reference_fasta: impl AsRef<Path>,
    contig_name_format: ContigNameFormat,
) -> Result<HashMap<String, Vec<u8>>, Error> {
    let fr = faidx::Reader::from_path(reference_fasta)?;

    let seq_names = fr.seq_names()?;
    let mut reference_genome_map = HashMap::with_capacity(seq_names.len());
    for seq_name in seq_names.iter() {
        reference_genome_map.insert(
            contig_name_format
                .contig_name_for_ref(&seq_name)
                .into_owned(),
            fr.fetch_seq(seq_name, 0, fr.fetch_seq_len(seq_name) as usize)?,
        );
    }

    Ok(reference_genome_map)
}

#[derive(Clone)]
pub struct GCCounter {
    mapq: u8,
    pub start_len: i64,
    pub end_len: i64,
    reference_genome_map: HashMap<String, Vec<u8>>,
    lag: usize,
    bam_path: String,
    contig_name_format: ContigNameFormat,
}

impl GCCounter {
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

    pub fn new(
        mapq: u8,
        start_len: i64,
        end_len: i64,
        reference_fasta: impl AsRef<Path>,
        lag: usize,
        bam_path: String,
    ) -> Result<Self, Error> {
        // load bam file header and see contig name starts with 'chr'
        let ir = IndexedReader::from_path(&bam_path)?;
        let header = ir.header();
        let contig_name_format = if header.target_names().iter().any(|t| t.starts_with(b"chr")) {
            ContigNameFormat::WithChr
        } else {
            ContigNameFormat::WithoutChr
        };

        event!(Level::DEBUG, "contig_name_format={contig_name_format:?}");

        let reference_genome_map = make_ref_seq_map(reference_fasta, contig_name_format)?;

        Ok(Self {
            mapq,
            start_len,
            end_len,
            reference_genome_map,
            lag,
            bam_path,
            contig_name_format,
        })
    }

    thread_local! {
        static BAM_MAP:RefCell<HashMap<String, IndexedReader>> = RefCell::new(
            HashMap::new()
        )

    }

    pub fn res_arr_shape(start_len: usize, end_len: usize) -> (usize, usize) {
        (end_len - start_len as usize + 1, 101)
    }

    pub fn make_res_arr(start_len: usize, end_len: usize) -> Array2<u64> {
        Array2::<u64>::zeros(Self::res_arr_shape(start_len, end_len))
    }

    pub fn count_gc(&self, fetch_region: &(String, i64, i64)) -> Result<Array2<u64>, Error> {
        let mut res_arr = Self::make_res_arr(self.start_len as usize, self.end_len as usize);

        Self::BAM_MAP.with_borrow_mut(|bam_map| {
            let ir = match bam_map.get_mut(&self.bam_path) {
                Some(ir) => ir,
                None => {
                    let ir = IndexedReader::from_path(&self.bam_path)?;
                    bam_map.insert(self.bam_path.clone(), ir);
                    bam_map.get_mut(&self.bam_path).unwrap()
                }
            };

            self._count_gc(fetch_region, ir, &mut res_arr)?;

            Ok::<_, Error>(())
        })?;

        Ok(res_arr)
    }

    #[inline]
    fn should_record_be_filtered_out(record: &Record) -> bool {
        record.is_unmapped()
            || record.is_mate_unmapped()
            || record.is_duplicate()
            || record.is_supplementary()
    }

    #[inline]
    fn ref_start_and_end_pos_for_gc(
        record: &Record,
        insert_size: i64,
        lag: usize,
    ) -> (usize, usize) {
        let ref_start = record.reference_start() as usize + lag;
        let ref_end = (record.reference_start() + insert_size) as usize - lag;

        (ref_start, ref_end)
    }

    #[inline]
    fn is_insert_size_not_target(insert_size: i64, start_len: i64, end_len: i64) -> bool {
        insert_size < start_len || insert_size > end_len
    }

    fn _count_gc(
        &self,
        fetch_region: &(String, i64, i64),
        ir: &mut IndexedReader,
        res_arr: &mut Array2<u64>,
    ) -> Result<(), Error> {
        let contig_name = &*self
            .contig_name_format
            .contig_name_for_ref(fetch_region.0.as_str());
        let contig_ref_seq = match self.reference_genome_map.get(contig_name) {
            Some(seq) => seq.as_slice(),
            None => Err(anyhow!("No contig name '{}' in ref seq map", contig_name))?,
        };

        match ir.fetch((contig_name, fetch_region.1, fetch_region.2)) {
            Ok(_) => {}
            Err(err) => Err(Error::from(err).context(format!(
                "{:?}",
                (fetch_region.0.as_str(), fetch_region.1, fetch_region.2,)
            )))?,
        };
        let mut record = Record::new();
        while let Some(rr) = ir.read(&mut record) {
            rr?;

            if Self::should_record_be_filtered_out(&record) {
                continue;
            }

            if record.mapq() < self.mapq || record.tid() != record.mtid() {
                continue;
            }

            let insert_size = record.insert_size();
            if Self::is_insert_size_not_target(insert_size, self.start_len, self.end_len) {
                continue;
            }

            let (ref_start, ref_end) =
                Self::ref_start_and_end_pos_for_gc(&record, insert_size, self.lag);

            if ref_end <= ref_start {
                // if insert size is very low. length <= 2*lag.
                continue;
            }

            let ref_seq = contig_ref_seq.get(ref_start..ref_end).ok_or_else(|| {
                anyhow!(
                    "Failed to get slice with index: {ref_start}..{ref_end}. contig={}",
                    contig_name
                )
            })?;

            if ref_seq.len() == 0 {
                continue;
            }

            let mut gc_cnt = 0;
            let mut at_cnt = 0;
            for &b in ref_seq {
                match Self::GC_TABLE[b as usize] {
                    1 => gc_cnt += 1,
                    2 => at_cnt += 1,
                    _ => {}
                }
            }

            let total_cnt = gc_cnt + at_cnt;
            let valid_percent = total_cnt as f64 / ref_seq.len() as f64;

            if valid_percent < 0.9 {
                continue;
            }

            // let gc_content = (gc_cnt as f64 / total_cnt as f64 * 100.0).round() as usize;
            let gc_content = integer_bankers_round(gc_cnt, total_cnt);

            match res_arr.get_mut(((insert_size - self.start_len) as usize, gc_content)) {
                Some(elem) => elem.add_assign(1),
                None => Err(anyhow!(
                    "Failed to get elem with index:({},{}) from result array.",
                    (insert_size - self.start_len),
                    gc_content
                ))?,
            }

            // println!("processed: {} {} {} {}", str::from_utf8(record.qname())?, gc_content, gc_cnt, total_cnt);
        }

        Ok(())
    }
}

#[inline]
fn python_round_logic(gc_cnt: i32, total_cnt: i32) -> usize {
    let ratio = gc_cnt as f64 / total_cnt as f64;
    let val_x100 = ratio * 100.0;

    let rounded = val_x100.round();
    let diff = rounded - val_x100;

    // Check if exactly halfway (within float epsilon)
    if (diff.abs() - 0.5).abs() < 1e-8 {
        // Tie-breaking: Round Half to Even
        if (rounded as i64) % 2 != 0 {
            // If rounded value is Odd (e.g. 13), force it to Even (12 or 14)
            // If diff > 0 (12.5 -> 13.0), we subtract 1 to get 12
            if diff > 0.0 {
                (rounded - 1.0) as usize
            } else {
                (rounded + 1.0) as usize
            }
        } else {
            // If rounded value is Even (e.g. 12), keep it
            rounded as usize
        }
    } else {
        rounded as usize
    }
}

#[inline]
fn integer_bankers_round(gc_cnt: i32, total_cnt: i32) -> usize {
    let numerator = 100 * gc_cnt;
    let denominator = total_cnt;

    let quotient = numerator / denominator;
    let remainder = numerator % denominator;
    let double_remainder = 2 * remainder;

    if double_remainder < denominator {
        quotient as usize
    } else if double_remainder > denominator {
        (quotient + 1) as usize
    } else {
        // Exact Tie (0.5) -> Round to Even
        if quotient % 2 == 0 {
            quotient as usize
        } else {
            (quotient + 1) as usize
        }
    }
}

fn make_correction_weight_arr(
    correction_weights_csv: impl AsRef<Path>,
    start_len: usize,
    end_len: usize,
) -> Result<Array2<f64>, Error> {
    let mut buf = String::new();
    let mut br = BufReader::new(File::open(correction_weights_csv)?);

    // parse first line
    if let Ok(true) = br.read_line(&mut buf).map(|n| n > 0) {
        let col_count = buf.trim().split(",").count();
        if col_count != 101 {
            Err(anyhow!("Expected 101 cols, but {col_count}"))?
        }
    }

    let n_rows = end_len - start_len + 1;
    let mut weights = Vec::with_capacity(101 * n_rows);
    let mut l_iter = start_len..(end_len + 1);

    while let Ok(true) = {
        buf.clear();
        br.read_line(&mut buf).map(|n| n > 0)
    } {
        l_iter.next().ok_or_else(|| {
            anyhow!(
                "Input file has extra lines exceeding the expected: {}.",
                end_len - start_len + 1
            )
        })?;

        for v in buf.trim().split(",") {
            weights.push({
                match v.parse::<f64>() {
                    Ok(f) => f,
                    Err(err) => Err(Error::from(err).context(v.to_string()))?,
                }
            });
        }
    }

    if l_iter.next().is_some() {
        Err(anyhow!(
            "Input file is too short. It does not cover the full range from {} to {}.",
            start_len,
            end_len
        ))?
    }

    Ok(Array2::from_shape_vec((n_rows, 101), weights)?)
}

/// Make an array of fragment sizes, and an array of GC-fragment size correction weights.
pub struct WeightedFragmentCollector {
    correction_weights_arr: ndarray::Array2<f64>,
    start_len: usize,
    end_len: usize,
    reference_seq_map: HashMap<String, Vec<u8>>,
    threads: usize,
    lag: usize,
}

impl WeightedFragmentCollector {
    const NA_WEIGHT: f64 = -1.0;
    const DEFAULT_CONTIG_NAME_FORMAT: ContigNameFormat = ContigNameFormat::WithChr;

    /// `start_len` and `end_len` are inclusive.
    pub fn new(
        correction_weights_csv: impl AsRef<Path>,
        start_len: usize,
        end_len: usize,
        lag: usize,
        reference_fasta: impl AsRef<Path>,
        threads: usize,
    ) -> Result<Self, Error> {
        event!(Level::WARN, "Filter logic is not implemented yet.");

        let reference_seq_map =
            make_ref_seq_map(reference_fasta, Self::DEFAULT_CONTIG_NAME_FORMAT)?;

        // build correction_weights_arr
        let correction_weights_arr =
            make_correction_weight_arr(correction_weights_csv, start_len, end_len)?;

        Ok(Self {
            correction_weights_arr,
            start_len,
            end_len,
            reference_seq_map,
            threads,
            lag,
        })
    }

    pub fn make_fs_and_correction_weight_arrs(
        &self,
        bam_path: impl AsRef<Path> + Sync + Send,
        region_infos: &[((&str, i64), &HashSet<&[u8]>)],
    ) -> Result<Vec<[(Array1<i32>, Array1<f64>); 2]>, Error> {
        let tp = ThreadPoolBuilder::new().num_threads(self.threads).build()?;

        let chunk_size = (region_infos.len() / self.threads).min(64).max(1);

        let ir_init = || IndexedReader::from_path(&bam_path);
        let res_vec = tp.install(|| {
            let r = region_infos
                .par_iter()
                .chunks(chunk_size)
                .map_init(ir_init, |ir_state, chunk| {
                    let ir = match ir_state {
                        Ok(v) => v,
                        Err(err) => Err(anyhow!(
                            "Failed to load bam: {err:?} {}",
                            bam_path.as_ref().display()
                        ))?,
                    };

                    let mut worker = WeightedFragmentCollectorWorker::default();

                    let mut res = Vec::with_capacity(chunk.len());

                    for &(region, separate_read_names) in chunk {
                        res.push(
                            worker.make_fs_and_correction_weight_arrs_for_region(
                                ir,
                                &self.correction_weights_arr,
                                region,
                                separate_read_names,
                                self.reference_seq_map
                                    .get(
                                        &*Self::DEFAULT_CONTIG_NAME_FORMAT
                                            .contig_name_for_ref(region.0),
                                    )
                                    .ok_or_else(|| {
                                        anyhow!(
                                            "Contig key not found in reference sequence map: {}",
                                            region.0
                                        )
                                    })?,
                                self.start_len,
                                self.end_len,
                                self.lag,
                            )?,
                        );
                    }

                    Ok::<_, Error>(res)
                })
                .collect::<Result<Vec<_>, Error>>()?
                .into_iter()
                .flatten()
                .collect::<Vec<_>>();

            Ok::<_, Error>(r)
        })?;

        Ok(res_vec)
    }
}

#[derive(Default)]
struct WeightedFragmentCollectorWorker {
    read_name_buf: Vec<u8>,
    group_a_fs_buf: Vec<i32>,
    group_a_weight_buf: Vec<f64>,
    group_b_fs_buf: Vec<i32>,
    group_b_weight_buf: Vec<f64>,
}

impl WeightedFragmentCollectorWorker {
    fn make_fs_and_correction_weight_arrs_for_region(
        &mut self,
        bam_reader: &mut IndexedReader,
        correction_weights_arr: &ndarray::Array2<f64>,
        region: (&str, i64),
        group_a_read_names: &HashSet<&[u8]>,
        contig_ref_seq: &[u8],
        start_len: usize,
        end_len: usize,
        lag: usize,
    ) -> Result<[(Array1<i32>, Array1<f64>); 2], Error> {
        let contig = region.0;
        let pos_0b = region.1 - 1;
        bam_reader.fetch((contig, pos_0b, pos_0b + 1))?;

        let read_name_buf = &mut self.read_name_buf;
        let group_a_fs_buf = &mut self.group_a_fs_buf;
        let group_a_weight_buf = &mut self.group_a_weight_buf;
        let group_b_fs_buf = &mut self.group_b_fs_buf;
        let group_b_weight_buf = &mut self.group_b_weight_buf;

        for pr in bam_reader.pileup_with_option(PileupOption {
            max_depth: i32::MAX,
            ignore_overlaps: true,
        }) {
            let plp = pr?;

            if (plp.pos() as i64) < pos_0b {
                continue;
            }

            for algmt in plp.alignments() {
                let record = algmt.record();

                read_name_buf.clear();
                read_name_buf.extend(record.qname());
                let suffix = if record.is_first_in_template() {
                    b"/1"
                } else {
                    b"/2"
                };
                read_name_buf.extend(suffix);

                let push_fs_and_weight = |fs_buf: &mut Vec<i32>, weight_buf: &mut Vec<f64>| {
                    let insert_size_abs = record.insert_size().abs();

                    let mut weight = WeightedFragmentCollector::NA_WEIGHT;

                    'weight_calc: {
                        let (ref_start, ref_end) =
                            GCCounter::ref_start_and_end_pos_for_gc(&record, insert_size_abs, lag);

                        if ref_end <= ref_start {
                            // if insert size is very low. length <= 2*lag.
                            break 'weight_calc;
                        }

                        let ref_seq = contig_ref_seq.get(ref_start..ref_end).ok_or_else(|| {
                            anyhow!(
                                "Failed to get slice with index: {ref_start}..{ref_end}. contig={}",
                                region.0
                            )
                        })?;

                        if ref_seq.len() == 0 {
                            break 'weight_calc;
                        }

                        let mut gc_cnt = 0;
                        let mut at_cnt = 0;
                        for &b in ref_seq {
                            match GCCounter::GC_TABLE[b as usize] {
                                1 => gc_cnt += 1,
                                2 => at_cnt += 1,
                                _ => {}
                            }
                        }

                        let total_cnt = gc_cnt + at_cnt;
                        let valid_percent = total_cnt as f64 / ref_seq.len() as f64;

                        if valid_percent < 0.9 {
                            break 'weight_calc;
                        }

                        // let gc_content = (gc_cnt as f64 / total_cnt as f64 * 100.0).round() as usize;
                        let gc_content = integer_bankers_round(gc_cnt, total_cnt);

                        let row_index = match (insert_size_abs as usize).checked_sub(start_len) {
                            Some(v) => v,
                            None => break 'weight_calc,
                        };

                        match correction_weights_arr.get((row_index, gc_content)).copied() {
                            Some(v) => weight = v,
                            None => Err(anyhow!(
                                "Failed to index weight with index: {:?}",
                                (row_index, gc_content)
                            ))?,
                        }
                    }

                    fs_buf.push(insert_size_abs as i32);
                    weight_buf.push(weight);

                    Ok::<_, Error>(())
                };

                if group_a_read_names.contains(read_name_buf.as_slice()) {
                    // always include this read. we won't check the rest filters.
                    push_fs_and_weight(group_a_fs_buf, group_a_weight_buf)?;
                } else {
                    // TODO: The same record filter logic as variant caller is needed.
                    if Self::should_record_be_filtered_out(&record) {
                        continue;
                    }
                    push_fs_and_weight(group_b_fs_buf, group_b_weight_buf)?;
                };
            }

            break;
        }

        Ok([
            (
                Array1::from_iter(group_a_fs_buf.drain(..)),
                Array1::from_iter(group_a_weight_buf.drain(..)),
            ),
            (
                Array1::from_iter(group_b_fs_buf.drain(..)),
                Array1::from_iter(group_b_weight_buf.drain(..)),
            ),
        ])
    }

    fn should_record_be_filtered_out(record: &Record) -> bool {
        false // TODO
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integer_bankers_round() {
        // --- TIE BREAKING (The .5 cases) ---

        // 1/8 = 12.5% -> Round to Even (12)
        // Python: round(12.5) -> 12
        assert_eq!(integer_bankers_round(1, 8), 12, "Failed 12.5% -> 12");

        // 3/8 = 37.5% -> Round to Even (38)
        // Python: round(37.5) -> 38
        assert_eq!(integer_bankers_round(3, 8), 38, "Failed 37.5% -> 38");

        // 5/8 = 62.5% -> Round to Even (62)
        // Python: round(62.5) -> 62
        assert_eq!(integer_bankers_round(5, 8), 62, "Failed 62.5% -> 62");

        // 7/8 = 87.5% -> Round to Even (88)
        // Python: round(87.5) -> 88
        assert_eq!(integer_bankers_round(7, 8), 88, "Failed 87.5% -> 88");

        // 1/40 = 2.5% -> Round to Even (2)
        // Python: round(2.5) -> 2
        assert_eq!(integer_bankers_round(1, 40), 2, "Failed 2.5% -> 2");

        // 3/40 = 7.5% -> Round to Even (8)
        // Python: round(7.5) -> 8
        assert_eq!(integer_bankers_round(3, 40), 8, "Failed 7.5% -> 8");

        // 1/200 = 0.5% -> Round to Even (0)
        // Python: round(0.5) -> 0
        assert_eq!(integer_bankers_round(1, 200), 0, "Failed 0.5% -> 0");

        // 3/200 = 1.5% -> Round to Even (2)
        // Python: round(1.5) -> 2
        assert_eq!(integer_bankers_round(3, 200), 2, "Failed 1.5% -> 2");

        // --- STANDARD ROUNDING ---

        // 1/3 = 33.33% -> 33
        assert_eq!(integer_bankers_round(1, 3), 33, "Failed 33.3%");

        // 2/3 = 66.66% -> 67
        assert_eq!(integer_bankers_round(2, 3), 67, "Failed 66.6%");

        // 1/2 = 50% (Exact)
        assert_eq!(integer_bankers_round(1, 2), 50, "Failed 50%");
    }

    #[test]
    fn test_python_rounding_parity() {
        // Case 1: 1/8 = 12.5%
        // Python: round(12.5) -> 12 (Round down to even)
        // Rust Default: 13 (Round up)
        assert_eq!(python_round_logic(1, 8), 12, "Failed 12.5% -> 12");

        // Case 2: 3/8 = 37.5%
        // Python: round(37.5) -> 38 (Round up to even)
        // Rust Default: 38 (Round up)
        assert_eq!(python_round_logic(3, 8), 38, "Failed 37.5% -> 38");

        // Case 3: 5/8 = 62.5%
        // Python: round(62.5) -> 62 (Round down to even)
        // Rust Default: 63 (Round up)
        assert_eq!(python_round_logic(5, 8), 62, "Failed 62.5% -> 62");

        // Case 4: 7/8 = 87.5%
        // Python: round(87.5) -> 88 (Round up to even)
        // Rust Default: 88 (Round up)
        assert_eq!(python_round_logic(7, 8), 88, "Failed 87.5% -> 88");

        // Case 5: 1/2 = 50.0% (Exact integer)
        assert_eq!(python_round_logic(1, 2), 50, "Failed 50%");

        // Case 6: 1/3 = 33.333...% (Standard rounding)
        assert_eq!(python_round_logic(1, 3), 33, "Failed 33.3%");

        // Case 7: 2/3 = 66.666...% (Standard rounding up)
        assert_eq!(python_round_logic(2, 3), 67, "Failed 66.6%");

        // Case 8: 0.5% (Very small) -> 0.005 * 100 = 0.5 -> rounds to 0 (even)
        // Constructing 0.5% is hard with integers, let's use 1/200 = 0.005 -> 0.5%
        assert_eq!(python_round_logic(1, 200), 0, "Failed 0.5% -> 0");

        // Case 9: 1.5% -> 3/200 = 0.015 -> 1.5% -> rounds to 2 (even)
        assert_eq!(
            python_round_logic(3, 200),
            2,
            "Faile
        d 1.5% -> 2"
        );
    }

    #[test]
    fn test_make_correction_weight_arr() -> Result<(), Box<dyn std::error::Error>> {
        let correction_weights_csv = "../../Sample_Output/Correction_Factors/test.csv";
        let start_len = 51;
        let end_len = 55;

        let arr = make_correction_weight_arr(correction_weights_csv, start_len, end_len)?;

        eprintln!("{arr:?}");

        Ok(())
    }
}
