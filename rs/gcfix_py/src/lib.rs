use std::{collections::HashSet, path::Path};

use gcfix::core::WeightedFragmentCollector;
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass(name = "WeightedFragmentCollector")]
pub(crate) struct WeightedFragmentCollectorPy {
    inner: WeightedFragmentCollector,
}

#[pymethods]
impl WeightedFragmentCollectorPy {
    #[new]
    fn new(
        correction_weights_csv: &str,
        start_len: usize,
        end_len: usize,
        lag: usize,
        reference_fasta: &str,
        threads: usize,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: WeightedFragmentCollector::new(
                correction_weights_csv,
                start_len,
                end_len,
                lag,
                reference_fasta,
                threads,
            )?,
        })
    }

    fn get_fs_and_weights<'py>(
        &self,
        py: Python<'py>,
        bam_path: &str,
        regions: Vec<(String, i64)>,
        group_a_read_names: Vec<Vec<String>>,
    ) -> PyResult<Vec<[(Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<f64>>); 2]>> {
        if regions.len() != group_a_read_names.len() {
            Err(PyValueError::new_err(format!(
                "Lengths of `regions` and `group_a_read_names` differ. {} != {}",
                regions.len(),
                group_a_read_names.len()
            )))?
        }

        let group_a_read_names_sets = group_a_read_names
            .iter()
            .map(|e| HashSet::from_iter(e.iter().map(String::as_bytes)))
            .collect::<Vec<_>>();

        let region_infos = regions
            .iter()
            .zip(group_a_read_names_sets.iter())
            .map(|(a, b)| ((a.0.as_str(), a.1), b))
            .collect::<Vec<_>>();

        let res = self
            .inner
            .make_fs_and_correction_weight_arrs(bam_path, &region_infos)?
            .into_iter()
            .map(|e| {
                let [(a0, a1), (b0, b1)] = e;
                [
                    (a0.into_pyarray(py), a1.into_pyarray(py)),
                    (b0.into_pyarray(py), b1.into_pyarray(py)),
                ]
            })
            .collect::<Vec<_>>();

        Ok(res)
    }
}

#[pymodule(name = "gcfix")]
mod core {
    #[pymodule_export]
    use super::WeightedFragmentCollectorPy;
}

#[cfg(test)]
mod tests {
    use crackle_kit::{
        tracing::{Level, event, level_filters::LevelFilter},
        tracing_kit::{setup_logging_stderr_only, setup_logging_stderr_only_debug},
    };
    use pyo3::{Python, ffi::c_str};

    use super::*;

    #[test]
    fn test_wfc_py() -> Result<(), Box<dyn std::error::Error>> {
        setup_logging_stderr_only_debug(LevelFilter::DEBUG)?;
        event!(Level::DEBUG, "Log enabled");
        pyo3::append_to_inittab!(core);

        Python::attach(|py| {
            // Python::run(py, c_str!(include_str!("../test_wfc_py.py")), None, None)
            Python::run(
                py,
                c_str!(
                    r##"
import pprint
import gcfix
from gcfix import WeightedFragmentCollector

import os
from pathlib import Path
print(Path.cwd().resolve())
wfc = WeightedFragmentCollector(
    correction_weights_csv="../../Sample_Output/Correction_Factors/sample1.csv",
    start_len=51,
    end_len=400,
    lag=10,
    reference_fasta="/home/eck/workspace/common_resources/hg38.fa",
    threads=8,
)

r = wfc.get_fs_and_weights(
    bam_path="../../Input_Bam/sample1.bam",
    regions=[
        ("1", 1139770),
        ("1", 1453962),
    ],
    group_a_read_names=[
        [
            "K00250:204:HNTGYBBXX:1:1203:30086:30626/1",
            # "K00250:204:HNTGYBBXX:1:1203:30086:30626/2",
        ],
        [
            "K00250:204:HNTGYBBXX:4:2223:22019:36798/1",
        ]
    ]
)

pprint.pprint(r)

                "##
                ),
                None,
                None,
            )
        })?;

        Ok(())
    }
}
