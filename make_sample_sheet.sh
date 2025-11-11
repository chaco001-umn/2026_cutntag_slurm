#!/bin/bash

mkdir -p /projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/
mkdir -p /projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/metadata
outfile=/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/metadata/samplesheet.csv

fastqdir=/projects/standard/swernimo/shared/Gearhart_Project_002_raw
R1s=( $(find "$fastqdir" -name "*R1_001.fastq.gz") )
rm -f $outfile

echo "sample    treatment:
SAW002	CT29_STB_040424_H3k9ac
SAW004	CT29_SR2_H3k9ac
SAW005	CT29_ST1_H3k9ac
SAW007	CT29_SR_DMSO_control_H3k9ac
SAW008	CT29_SR_DMSO_control_H3k9ac
SAW009	CT29_SR_UK_H3k9ac
SAW010	CT29_SR_UK_H3k9ac
SAW011	CT29_ST_DMSO_control_H3k9ac
SAW012	CT29_ST_DMSO_control_H3k9ac
SAW013	CT29_ST_UK_H3k9ac"

for R1 in ${R1s[@]}; do
    case $R1 in
        *SAW002*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_STB,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW004*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_SR2,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW005*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_ST1,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW007*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_SR_DMSO_control,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW008*)
            rep=2
            R2=${R1/R1_001/R2_001}
            echo "CT29_SR_DMSO_control,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW009*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_SR_UK_H3k9ac,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW010*)
            rep=2
            R2=${R1/R1_001/R2_001}
            echo "CT29_SR_UK_H3k9ac,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW011*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_ST_DMSO_contro,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW012*)
            rep=2
            R2=${R1/R1_001/R2_001}
            echo "CT29_ST_DMSO_contro,$rep,$R1,$R2" >> $outfile
        ;;
        *SAW013*)
            rep=1
            R2=${R1/R1_001/R2_001}
            echo "CT29_ST_UK,$rep,$R1,$R2" >> $outfile
        ;;
    esac
done

head -1000 $outfile

exit 0