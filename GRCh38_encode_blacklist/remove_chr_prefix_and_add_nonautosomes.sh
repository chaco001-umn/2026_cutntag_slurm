#!/bin/bash

input_bed=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/GRCh38_encode_blacklist/ENCFF356LFX.bed
output_bed=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/GRCh38_encode_blacklist/blacklist.bed

sed 's/^chr//' ${input_bed} > ${output_bed}

ref=/common/bioref/ensembl/main/Homo_sapiens-113/GRCh38.p14/seq/genome.fa

cut -f1,2 ${ref}.fai | \
    tail -n+25 | \
    awk 'BEGIN{OFS="\t"} {print $1, 1, $2}' >> ${output_bed}

