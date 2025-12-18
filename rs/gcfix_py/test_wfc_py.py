import gcfix
from gcfix import WeightedFragmentCollector

import os
print(os.curdir())
WeightedFragmentCollector(
    correction_weights_csv="../../Sample_Output/Correction_Factors/sample1.csv",
    start_len=51,
    end_len=400,
    lag=10,
    reference_fasta="/home/eck/workspace/common_resources/hg38.fa",
    threads=8,
)
