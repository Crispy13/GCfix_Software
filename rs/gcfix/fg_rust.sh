set -e

export LD_LIBRARY_PATH="/home/eck/software/miniconda3/envs/gcfix/lib/"

input_bam="../NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
out_npy="$(basename $input_bam).npy"

command time -v \
cargo flamegraph -o fg.svg --profile debug-release -- \
    $input_bam \
    $out_npy \
    --mapq 30 \
    -s 51 \
    -e 400 \
    -r /home/eck/workspace/common_resources/hg38.fa \
    -b ../hg38/GC_tagging_bin_locations.csv \
    -t 8 \
    --lag 10

# ../Input_Bam/sample1.bam \