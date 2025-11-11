#!/bin/bash

apptainer_dir=/projects/standard/swernimo/shared/RI/software/apptainers

bind_paths=${1}

fastp(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/fastp.sif fastp "$@"
}
export -f fastp

bowtie2(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/bowtie2.sif bowtie2 "$@"
}
export -f bowtie2

picard(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/picard.sif picard "$@"
}
export -f picard

samtools(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/samtools.sif samtools "$@"
}
export -f samtools

macs3(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/macs3.sif macs3 "$@"
}
export -f macs3

Genrich(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/genrich.sif Genrich "$@"
}
export -f Genrich

multiqc(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/multiqc.sif multiqc "$@"
}
export -f multiqc

bedtools(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/bedtools.sif bedtools "$@"
}
export -f bedtools

seacr(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/seacr.sif /opt/conda/envs/bio/bin/SEACR_1.3.sh "$@"
}
export -f seacr

gopeaks(){
    apptainer exec -B ${bind_paths} ${apptainer_dir}/gopeaks.sif gopeaks "$@"
}
export -f gopeaks