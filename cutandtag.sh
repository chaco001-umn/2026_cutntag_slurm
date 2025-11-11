#!/bin/bash

set -euo pipefail

export functions_file=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/bash_functions.sh
export single_sample_script=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/single_sample.sbatch
export multi_sample_script=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/multi_sample.sbatch
# this blacklist was made by taking encode's list, then adding MT and all non-autosomes in their full length
export GRCh38_blacklist_bed=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/GRCh38_encode_blacklist/blacklist.bed
source ${functions_file}

#### Project-specific:
export outdir=/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/
export samplesheet=${outdir}/metadata/samplesheet.csv
export bowtie2_index=/common/bioref/ensembl/main/Homo_sapiens-113/GRCh38.p14/bowtie2/genome
export genome_fai=/common/bioref/ensembl/main/Homo_sapiens-113/GRCh38.p14/seq/genome.fa.fai

n_samples=$(get_n_samples $samplesheet)
fqpath=$(get_fq_path $samplesheet)

echo "$n_samples samples found in sample sheet"
echo "fastq files are here: ${fqpath}"
arraymax=$(expr $n_samples - 1)

export bind_paths="${fqpath}:${fqpath},${outdir}:${outdir},/scratch.global/chaco001:/tmp,/common/bioref:/common/bioref"
export apptainers_file=/projects/standard/swernimo/shared/RI/2025_cutntag/pipeline/apptainer_paths.sh
source ${apptainers_file} "${bind_paths}"


make_output_dirs $outdir

# check for the single sample done files
n_completed_samples=$(count_single_sample_done_files $outdir)
echo "$n_completed_samples have already completed"

if [ ! "$n_completed_samples" = "$n_samples" ]; then
        echo "submitting single sample job array and multi-sample job..."
    single_sample_job=$(sbatch \
        --array=0-$arraymax \
        --mail-type="begin,end,fail" \
        --mail-user "$USER@umn.edu" \
        -J "singlesample" \
        -o "${outdir}/slurmlogs/singlesample_%A_%a.out" \
        -e "${outdir}/slurmlogs/singlesample_%A_%a.err" \
        ${single_sample_script})
    single_sample_jobid=$(echo "$single_sample_job" | awk '{print $4}')

    multi_sample_job=$(sbatch \
        --mail-type="begin,end,fail" \
        --mail-user "$USER@umn.edu" \
        -J "multisample" \
        -o "${outdir}/slurmlogs/multisample_%A_%a.out" \
        -e "${outdir}/slurmlogs/multisample_%A_%a.err" \
        --dependency=afterok:${single_sample_jobid} \
        --kill-on-invalid-dep=yes \
        ${multi_sample_script})

else
    echo "all single sample jobs finished, only submitting multi-sample job..."
    multi_sample_job=$(sbatch \
        --mail-type="begin,end,fail" \
        --mail-user "$USER@umn.edu" \
        -J "multisample" \
        -o "${outdir}/slurmlogs/multisample_%A_%a.out" \
        -e "${outdir}/slurmlogs/multisample_%A_%a.err" \
        ${multi_sample_script})
fi
exit 0



