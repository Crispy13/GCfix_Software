use std::{cell::RefCell, collections::HashMap, ops::AddAssign, os::unix::thread, path::{Path, PathBuf}};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::locus::GenomeRegion,
    rust_htslib::bam::{IndexedReader, Read, Record, ext::BamRecordExtensions as _},
};
use ndarray::{Array2, ArrayBase, Dim, OwnedRepr};

pub struct GCCounter {
    mapq: u8,
    start_len: i64,
    end_len: i64,
    reference_genome_map: HashMap<Vec<u8>, Vec<u8>>,
    lag: usize,
    bam_path: String,
}

impl GCCounter {
    pub fn new(
        mapq: u8,
        start_len: i64,
        end_len: i64,
        reference_fasta: impl AsRef<Path>,
        lag: usize,
        bam_path: String,
    ) -> Self {
        Self {
            mapq,
            start_len,
            end_len,
            reference_genome_map,
            lag,
            bam_path,
        }
    }

    thread_local! {
        static BAM_MAP:RefCell<HashMap<String, IndexedReader>> = RefCell::new(
            HashMap::new()
        )

    }

    pub fn count_gc(&self, fetch_region: &(&[u8], i64, i64)) -> Result<(), Error> {
        let mut res_arr =
            Array2::<u64>::zeros((self.end_len as usize - self.start_len as usize + 1, 101));
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

        Ok(())
    }

    fn _count_gc(
        &self,
        fetch_region: &(&[u8], i64, i64),
        ir: &mut IndexedReader,
        res_arr: &mut Array2<u64>,
    ) -> Result<(), Error> {
        let contig_name = fetch_region.0;
        let contig_ref_seq = match self.reference_genome_map.get(contig_name) {
            Some(seq) => seq.as_slice(),
            None => Err(anyhow!(
                "No contig name '{}' in ref seq map",
                String::from_utf8_lossy(contig_name)
            ))?,
        };

        ir.fetch(*fetch_region)?;
        let mut record = Record::new();
        while let Some(rr) = ir.read(&mut record) {
            rr?;

            if record.is_unmapped()
                || record.is_mate_unmapped()
                || record.is_duplicate()
                || record.is_supplementary()
            {
                continue;
            }

            if record.mapq() < self.mapq || record.tid() != record.mtid() {
                continue;
            }

            let length = record.insert_size();
            if length < self.start_len || length > self.end_len {
                continue;
            }

            let ref_start = record.reference_start() as usize + self.lag;
            let ref_end = (record.reference_start() + length) as usize - self.lag;
            let ref_seq = contig_ref_seq.get(ref_start..ref_end).ok_or_else(|| {
                anyhow!(
                    "Failed to get slice with index: {ref_start}..{ref_end}. contig={}",
                    String::from_utf8_lossy(contig_name)
                )
            })?;

            if ref_seq.len() == 0 {
                continue;
            }

            let mut gc_cnt = 0;
            let mut at_cnt = 0;
            for b in ref_seq {
                match b.to_ascii_lowercase() {
                    b'g' | b'c' => {
                        gc_cnt += 1;
                    }
                    b'a' | b't' => {
                        at_cnt += 1;
                    }
                    _oth => continue,
                }
            }

            let total_cnt = gc_cnt + at_cnt;
            let valid_percent = total_cnt as f64 / ref_seq.len() as f64;

            if valid_percent < 0.9 {
                continue;
            }

            let gc_content = (gc_cnt as f64 / total_cnt as f64 * 100.0).round() as usize;

            match res_arr.get_mut(((length - self.start_len) as usize, gc_content)) {
                Some(elem) => elem.add_assign(1),
                None => Err(anyhow!(
                    "Failed to get elem with index:({},{}) from result array.",
                    (length - self.start_len),
                    gc_content
                ))?,
            }
        }

        Ok(())
    }
}
