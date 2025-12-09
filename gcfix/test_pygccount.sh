set -e

input_bam="../NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
out_npy="$(basename $input_bam).py.npy"

command time -v \
python ../single_sample_GC_count.py \
    $input_bam \
    $out_npy \
    30 \
    51 \
    400 \
    8 \
    /home/eck/workspace/common_resources/hg38.fa \
    ../hg38/GC_tagging_bin_locations.csv \

# ../Input_Bam/sample1.bam \