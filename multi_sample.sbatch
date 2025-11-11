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

# make the initial multiqc report
if [ ! -f ${outdir}/reports/single_sample_qc.html ]; then
    multiqc --outdir ${outdir}/reports \
        --filename single_sample_qc.html \
        --title "Single Sample read and mapping qc" \
        ${outdir}/qc_singlesample/fastp \
        ${outdir}/qc_singlesample/samtools \
        ${outdir}/qc_singlesample/picard
fi

# do additional QC that multiqc doesn't know how to handle

#####
# compile fraction-reads-in-peaks into single files
####
cd ${outdir}/macs3
samples=($(find . -maxdepth 1 -type f -name "*_peaks.narrowPeak" | sed 's|^\./||'))
echo -e "sample\tn_peaks\tn_reads_in_peaks\tn_mapped_reads\tFRiP" > $outdir/qc_multisample/macs3_FRiP.txt
for sample in ${samples[@]}; do
    tail -1 ${sample/_peaks.narrowPeak/_reads_in_peaks.tsv} | \
        awk '{$1=$1; print}' OFS="\t" | \
        awk '{printf "%s\t%.6f\n", $0, $3/$4}' >> $outdir/qc_multisample/macs3_FRiP.txt

done

cd ${outdir}/genrich
samples=($(find . -maxdepth 1 -type f -name "*_peaks.narrowPeak" | sed 's|^\./||'))
echo -e "sample\tn_peaks\tn_reads_in_peaks\tn_mapped_reads\tFRiP" > $outdir/qc_multisample/genrich_FRiP.txt
for sample in ${samples[@]}; do
    tail -1 ${sample/_peaks.narrowPeak/_reads_in_peaks.tsv} | \
        awk '{$1=$1; print}' OFS="\t" | \
        awk '{printf "%s\t%.6f\n", $0, $3/$4}' >> $outdir/qc_multisample/genrich_FRiP.txt

done
#####
# find out how much within-caller overlap there is among samples
#####
cd ${outdir}/macs3
samples=($(find . -maxdepth 1 -type f -name "*_peaks.narrowPeak" | sed 's|^\./||'))
n=${#samples[@]}
echo -e "sample1\tsample2\tsample1_peaks\tsample2_peaks\tintersection_bases\tunion_bases\tjaccard\tn_intersections\tn_50_intersections" > $outdir/qc_multisample/macs3_intersample_peak_overlap_metrics.txt
for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
        sample1="${samples[i]}"
        sample2="${samples[j]}"
        sample1_peaks=$(wc -l ${sample1} | cut -f1 -d" ")
        sample2_peaks=$(wc -l ${sample2} | cut -f1 -d" ")
        half_overlap=$(bedtools intersect -u -f 0.50 -r -a "$sample1" -b "$sample2" | wc -l)
        jaccard=$(bedtools jaccard -a "$sample1" -b "$sample2" | \
                tail -1 | awk '{$1=$1; print}' OFS="\t")
        
        echo -e "${sample1/_peaks.narrowPeak/}\t${sample2/_peaks.narrowPeak/}\t${sample1_peaks}\t${sample2_peaks}\t${jaccard}\t${half_overlap}" >> $outdir/qc_multisample/macs3_intersample_peak_overlap_metrics.txt
    done
done

# genrich
cd ${outdir}/genrich
samples=($(find . -maxdepth 1 -type f -name "*_peaks.narrowPeak" | sed 's|^\./||'))
n=${#samples[@]}
echo -e "sample1\tsample2\tsample1_peaks\tsample2_peaks\tintersection_bases\tunion_bases\tjaccard\tn_intersections\tn_50_intersections" > $outdir/qc_multisample/genrich_intersample_peak_overlap_metrics.txt
for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
        sample1="${samples[i]}"
        sample2="${samples[j]}"
        sample1_peaks=$(wc -l ${sample1} | cut -f1 -d" ")
        sample2_peaks=$(wc -l ${sample2} | cut -f1 -d" ")        
        half_overlap=$(bedtools intersect -u -f 0.50 -r -a "$sample1" -b "$sample2" | wc -l)
        jaccard=$(bedtools jaccard -a "$sample1" -b "$sample2" | \
                tail -1 | awk '{$1=$1; print}' OFS="\t")
        
        echo -e "${sample1/_peaks.narrowPeak/}\t${sample2/_peaks.narrowPeak/}\t${sample1_peaks}\t${sample2_peaks}\t${jaccard}\t${half_overlap}">> $outdir/qc_multisample/genrich_intersample_peak_overlap_metrics.txt
    done
done

# macs3 vs genrich, sample by sample (as opposed to intersample)
echo -e "sample\tn_macs3_peaks\tn_genrich_peaks\tintersection\tunion\tjaccard\tn_intersections\tn_50_intersections" > $outdir/qc_multisample/macs3_vs_genrich_intrasample_peak_overlap_metrics.txt

for sample in ${samples[@]}; do
    sample_peaks_macs3=$(wc -l < "${outdir}/macs3/${sample}")
    sample_peaks_genrich=$(wc -l < "${outdir}/genrich/${sample}")
    half_overlap=$(bedtools intersect -u -f 0.50 -r -a "${outdir}/macs3/$sample" -b "${outdir}/genrich/$sample" | wc -l)
    jaccard=$(bedtools jaccard -a "${outdir}/macs3/$sample" -b "${outdir}/genrich/$sample" | \
                tail -1 | awk '{$1=$1; print}' OFS="\t")
    echo -e "${sample/_peaks.narrowPeak/}\t${sample_peaks_macs3}\t${sample_peaks_genrich}\t${jaccard}\t${half_overlap}" >> $outdir/qc_multisample/macs3_vs_genrich_intrasample_peak_overlap_metrics.txt
done


### Now get consensuses across samples (within treatments), then union across treatments


### then feature counts


exit 0


