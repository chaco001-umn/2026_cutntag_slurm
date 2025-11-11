#!/bin/bash


count_single_sample_done_files(){
    local outdir=${1}
    n_files=$(find "$outdir/checkpoints/" -maxdepth 1 -type f -name "*singlesample*" | wc -l)
    echo "$n_files"
}

get_n_samples(){
    local samplesheet=${1}
    local lines
    lines=$(wc -l < "$samplesheet")
    echo "$lines"
}
export -f get_n_samples

get_fq_path(){
    local samplesheet=${1}
    local fqpath
    fqpath=$(head -1 "$samplesheet" | cut -f3 -d,)
    fqpath=$(dirname $fqpath)
    echo "$fqpath"
}
export -f get_fq_path

make_output_dirs(){
    local outdir=${1}
    mkdir -p $outdir/fastp
    mkdir -p $outdir/bowtie2
    mkdir -p $outdir/dedup
    mkdir -p $outdir/macs3
    mkdir -p $outdir/genrich
    mkdir -p $outdir/seacr
    mkdir -p $outdir/gopeaks
    mkdir -p $outdir/checkpoints
    mkdir -p $outdir/slurmlogs
    mkdir -p $outdir/qc_singlesample
    mkdir -p $outdir/qc_multisample
    mkdir -p $outdir/reports
}
export -f make_output_dirs

get_samplename(){
    local samplesheet=${1}
    local samplenum=${2}
    local fqpath
    sample=$(head -"$samplenum" "$samplesheet" | tail -1 | cut -f3 -d,)
    sample=$(basename $sample)
    sample=${sample/_R1_001.fastq.gz}
    echo "$sample"
}
export -f get_samplename

get_R1(){
    local samplesheet=${1}
    local samplenum=${2}
    local fqpath
    R1=$(head -"$samplenum" "$samplesheet" | tail -1 | cut -f3 -d,)
    R1=$(basename $R1)
    echo "$R1"
}
export -f get_R1

get_R2(){
    local samplesheet=${1}
    local samplenum=${2}
    local fqpath
    R2=$(head -"$samplenum" "$samplesheet" | tail -1 | cut -f4 -d,)
    R2=$(basename $R2)
    echo "$R2"
}
export -f get_R2