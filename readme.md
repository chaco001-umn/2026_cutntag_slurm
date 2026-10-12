This is a slurm-based pipeline to analyze CUT&TAG data.

One needs to make all of the apptainers in apptainer_defs, and put them into the 
directory specified in apptainers.sh (that can be changed if desired).

Then a sample sheet is needed that looks like this:

Please note the paths should actually be full paths, starting with "/"
```
CT29_STB,1,SAW002_S70_R1_001.fastq.gz,SAW002_S70_R2_001.fastq.gz
CT29_SR2,1,SAW004_S71_R1_001.fastq.gz,SAW004_S71_R2_001.fastq.gz
CT29_ST1,1,SAW005_S72_R1_001.fastq.gz,SAW005_S72_R2_001.fastq.gz
```

Then, adjust relevant variables (like paths) in cutntag.sh, and run it. Or submit
 the single_sample and multi_sample sbatch scripts in order.
 
The pipeline runs fastp to trim, aligns with bowtie2, checks some
mapping statistics and removes blacklisted regions, marks duplicates,
calls peaks with macs3, genrich, seacr, and gopeaks, in the single sample stage.

The multi-sample stage calculates correlations between samples within callers,
and between callers (macs3, genrich only, due to discussion with PI) within samples.
Then it makes intersection bed files if at least two samples, within a caller,
detected a peak. We decided to keep the intersections separate by caller. 

Finally featureCounts is run for each intersected peak list. 


To produce statistics and figures, see the files in "example_R_analysis." These
were specific to the exact project, but are likely modifyable. Start with the
"...analysis.R" which makes objects and saves them, then proceed to "...statistics.R"

The Rmd files generated reports for the specific project. 