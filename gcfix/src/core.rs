use std::{
    borrow::Cow,
    cell::RefCell,
    collections::HashMap,
    ops::AddAssign,
    os::unix::thread,
    path::{Path, PathBuf},
};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::locus::GenomeRegion,
    rust_htslib::{
        bam::{IndexedReader, Read, Record, ext::BamRecordExtensions as _},
        faidx,
    },
};
use ndarray::{Array, Array2, ArrayBase, Dim, OwnedRepr};

#[derive(Clone)]
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
        let contig_name_format = if header.target_names().iter().any(|t| t.starts_with(b"chr"))
        {
            ContigNameFormat::WithChr
        } else {
            ContigNameFormat::WithoutChr
        };

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

    fn _count_gc(
        &self,
        fetch_region: &(String, i64, i64),
        ir: &mut IndexedReader,
        res_arr: &mut Array2<u64>,
    ) -> Result<(), Error> {
        let contig_name = &*self.contig_name_format.contig_name_for_ref(fetch_region.0.as_str());
        let contig_ref_seq = match self
            .reference_genome_map
            .get(contig_name)
        {
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
                    contig_name
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
