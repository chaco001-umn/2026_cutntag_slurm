#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=chaco001@umn.edu 
#SBATCH --time=4:00:00         
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=24      
#SBATCH --mem=128G               
#SBATCH --tmp=300G
#SBATCH --partition=msismall

set -euo


source ${functions_file}
source ${apptainers_file} "${bind_paths}"

sample_num=$(expr $SLURM_ARRAY_TASK_ID + 1)
echo ${sample_num}

fqpath=$(get_fq_path $samplesheet)
samplename=$(get_samplename ${samplesheet} ${sample_num})
R1=$(get_R1 ${samplesheet} ${sample_num})
R2=$(get_R2 ${samplesheet} ${sample_num})

echo "working on ${samplename} in ${fqpath}"
echo "R1 = $R1"
echo "R2 = $R2"

if [ ! -f ${outdir}/checkpoints/fastp_${samplename}.done ]; then

    echo "running fastp..."
    
    mkdir -p ${outdir}/qc_singlesample/fastp
    
    fastp -i ${fqpath}/${R1} -I ${fqpath}/${R2} \
        -o  ${outdir}/fastp/${R1} -O ${outdir}/fastp/${R2} \
        --length_required 10 \
        --detect_adapter_for_pe \
        --json ${outdir}/qc_singlesample/fastp/${samplename}.json

    touch ${outdir}/checkpoints/fastp_${samplename}.done
else
    echo "fastp .done file found, skipping..."
fi

if [ ! -f ${outdir}/checkpoints/bowtie2_${samplename}.done ]; then
    echo "running bowtie2..."
    
    mkdir -p ${outdir}/qc_singlesample/samtools
    
    # note: bowtie2_index is exported from cutntag.sh
    bowtie2 --local --very-sensitive-local --no-mixed --no-discordant \
      --threads 24 \
      -x ${bowtie2_index} -1 ${outdir}/fastp/${R1} -2 ${outdir}/fastp/${R2} \
        --rg-id "${samplename}" \
        --rg "SM:${samplename}" \
        --rg "PL:ILLUMINA" \
        --rg "LB:${samplename}" \
      -S ${outdir}/bowtie2/${samplename}_allcontigs.sam
    samtools view -bS ${outdir}/bowtie2/${samplename}_allcontigs.sam | \
        samtools sort -o ${outdir}/bowtie2/${samplename}_allcontigs.bam
    samtools index ${outdir}/bowtie2/${samplename}_allcontigs.bam
    
    rm ${outdir}/bowtie2/${samplename}_allcontigs.sam
    
    samtools flagstat ${outdir}/bowtie2/${samplename}_allcontigs.bam > ${outdir}/qc_singlesample/samtools/${samplename}_flagstats.txt
    samtools stat ${outdir}/bowtie2/${samplename}_allcontigs.bam | grep "^SN" > ${outdir}/qc_singlesample/samtools/${samplename}_SNstats.txt
    samtools view -f 3 ${outdir}/bowtie2/${samplename}_allcontigs.bam | \
        awk '{print $9}' | \
        awk '$1 > 0 && $1 < 1000' | \
        sort | uniq -c > ${outdir}/qc_singlesample/samtools/${samplename}_fragmenthist.txt
    samtools idxstats ${outdir}/bowtie2/${samplename}_allcontigs.bam | head -25 > ${outdir}/qc_singlesample/samtools/${samplename}_idxstats.txt
        
    # toss blacklist, non-autosomes/sex chroms
    bedtools intersect -v \
        -abam ${outdir}/bowtie2/${samplename}_allcontigs.bam \
        -b ${GRCh38_blacklist_bed} > ${outdir}/bowtie2/${samplename}.bam
    samtools index ${outdir}/bowtie2/${samplename}.bam


    touch ${outdir}/checkpoints/bowtie2_${samplename}.done
else
    echo "bowtie2 .done file found, skipping..."
fi

if [ ! -f ${outdir}/checkpoints/dedup_${samplename}.done ]; then
    echo "running picard to deduplicate..."
    
    mkdir -p ${outdir}/qc_singlesample/picard
    
    if [ ! -f ${outdir}/dedup/${samplename}.bam ]; then
    
        picard MarkDuplicates I=${outdir}/bowtie2/${samplename}.bam \
            O=${outdir}/dedup/${samplename}.bam \
            M=${outdir}/qc_singlesample/picard/${samplename}_picardmetrics.txt \
            REMOVE_DUPLICATES=true
    
    fi
        
    picard CollectInsertSizeMetrics I=${outdir}/bowtie2/${samplename}.bam \
        O=${outdir}/qc_singlesample/picard/${samplename}_insertsizemetrics.txt \
        H=${outdir}/qc_singlesample/picard/${samplename}_insertsizehist.pdf 
    
    touch ${outdir}/checkpoints/dedup_${samplename}.done
else
    echo "picard dedup .done file found, skipping..."
fi
    

if [ ! -f ${outdir}/checkpoints/macs3_${samplename}.done ]; then
    echo "running macs3..."
    
    # --shift, --extsize not needed because pe
    # --nolambda to use global background, not local "With this flag on, MACS3 will use the background lambda as local lambda. This means MACS3 will not consider the local bias at peak candidate regions. It is particularly recommended while calling peaks without control sample."
    macs3 callpeak -t ${outdir}/dedup/${samplename}.bam \
        -f BAMPE \
        -g hs \
        --outdir ${outdir}/macs3 \
        -n ${samplename} \
        --nomodel \
        --max-gap 300 \
        --cutoff-analysis \
        --extsize 200 \
        --keep-dup all \
        --qvalue 0.1
        
    # original
    #macs3 callpeak -t ${outdir}/dedup/${samplename}.bam \
    #    -f BAMPE \
    #    -g hs \
    #    --outdir ${outdir}/macs3 \
    #    -n ${samplename} \
    #    --nomodel \
    #    --max-gap 300 \
    #    --qvalue 0.05
        
    # do QC 
    total_mapped=$(samtools view -F 4 ${outdir}/dedup/${samplename}.bam | wc -l)
    in_peaks=$(samtools view -F 4 -L ${outdir}/macs3/${samplename}_peaks.narrowPeak ${outdir}/dedup/${samplename}.bam | wc -l)
    n_peaks=$(wc -l ${outdir}/macs3/${samplename}_peaks.narrowPeak | cut -f1 -d" ")
    echo "samplename    n_peaks n_reads_in_peaks    n_mapped_reads" > ${outdir}/macs3/${samplename}_reads_in_peaks.tsv
    echo "${samplename} ${n_peaks} ${in_peaks} ${total_mapped}" >> ${outdir}/macs3/${samplename}_reads_in_peaks.tsv
    

    touch ${outdir}/checkpoints/macs3_${samplename}.done
else
    echo "macs3 .done file found, skipping..."
fi

if [ ! -f ${outdir}/checkpoints/genrich_${samplename}.done ]; then
    echo "running genrich..."
    
    # first must sort by query name
    if [ ! -f ${outdir}/bowtie2/${samplename}.namesort.bam ]; then
        samtools sort -n ${outdir}/dedup/${samplename}.bam -o ${outdir}/bowtie2/${samplename}.namesort.bam
    fi
    
    # note: we dont bother with duplicate removal because
    # we already removed them with picard
    Genrich -t ${outdir}/bowtie2/${samplename}.namesort.bam \
        -o ${outdir}/genrich/${samplename}_peaks.narrowPeak \
        -f ${outdir}/genrich/${samplename}.log \
        -j -y \
        -q 0.25
    
    # original    
    #    Genrich -t ${outdir}/bowtie2/${samplename}.namesort.bam \
    #    -o ${outdir}/genrich/${samplename}_peaks.narrowPeak \
    #    -f ${outdir}/genrich/${samplename}.log \
    #    -j -y -v -e chrM
        
    # do QC 
    total_mapped=$(samtools view -F 4 ${outdir}/dedup/${samplename}.bam | wc -l)
    in_peaks=$(samtools view -F 4 -L ${outdir}/genrich/${samplename}_peaks.narrowPeak ${outdir}/dedup/${samplename}.bam | wc -l)
    n_peaks=$(wc -l ${outdir}/genrich/${samplename}_peaks.narrowPeak | cut -f1 -d" ")
    echo "samplename    n_peaks n_reads_in_peaks    n_mapped_reads" > ${outdir}/genrich/${samplename}_reads_in_peaks.tsv
    echo "${samplename} ${n_peaks} ${in_peaks} ${total_mapped}" >> ${outdir}/genrich/${samplename}_reads_in_peaks.tsv
    
    touch ${outdir}/checkpoints/genrich_${samplename}.done
else
    echo "genrich .done file found, skipping..."
fi


if [ ! -f ${outdir}/checkpoints/seacr_${samplename}.done ]; then
    echo "running seacr..."
    
    # need to keep just mapped pairs otherwise the bed file has "." and bedtools genomecov doesn't know how to deal
    samtools view -f 3 -b ${outdir}/dedup/${samplename}.bam > ${outdir}/seacr/${samplename}.bam
    # directly from the github; conversion to fragments bedgraph

    bedtools bamtobed -bedpe -i ${outdir}/seacr/${samplename}.bam > ${outdir}/seacr/${samplename}.bed
    cut -f 1,2,6 ${outdir}/seacr/${samplename}.bed | sort -k1,1 -k2,2n -k3,3n > ${outdir}/seacr/${samplename}.fragments.bed
    bedtools genomecov -bg -i ${outdir}/seacr/${samplename}.fragments.bed -g ${genome_fai} > ${outdir}/seacr/${samplename}.fragments.bedgraph
    
    seacr ${outdir}/seacr/${samplename}.fragments.bedgraph \
        0.1 \
        norm \
        relaxed \
        ${outdir}/seacr/${samplename}
    
        
    touch ${outdir}/checkpoints/seacr_${samplename}.done
else
    echo "seacr .done file found, skipping..."
fi

if [ ! -f ${outdir}/checkpoints/gopeaks_${samplename}.done ]; then
    echo "running gopeaks..."

    gopeaks -b ${outdir}/dedup/${samplename}.bam \
        -p 0.1 \
        -o ${outdir}/gopeaks/${samplename}
    
        
    touch ${outdir}/checkpoints/gopeaks_${samplename}.done
else
    echo "gopeaks .done file found, skipping..."
fi




echo "all steps complete, creating  ${outdir}/checkpoints/singlesample_${samplename}.done"
touch ${outdir}/checkpoints/singlesample_${samplename}.done
exit 0