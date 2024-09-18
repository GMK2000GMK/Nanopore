#!/bin/bash

# TB_Workflow
# This workflow is designed to analyze nanopore data in Linux. 
# You are expected to provide a reference sequence and a fastq file. 
# The fastq files are classified using a Centrifuge database. 
# The results are then filtered for Mycobacterium read IDs. 
# The filtered Mycobacterium fastq files are mapped to your selected reference 
# and converted into a SAM file, then BAM file, and finally a sorted BAM file using Samtools. 
# The coverage at each position is obtained. 
# Once these processes are complete, the important files are run through a script 
# which outputs a CSV file with the calculated results from the run. 
# Finally, from the sorted BAM file, a Coverage plot and VCF file are produced.

# Variables
run=''
barcode=''
reference=''
database=''

# Create necessary directories
cd ~
mkdir -p ONT_ANALYSIS/{fastq,temp_files1,important_files,output}

# Assign the reference, barcode, and run name to create a combined single fastq file
echo "Barcode: $barcode"
echo "Run: $run"
files="/fastq_all/"
data=$(ls -m ${files} | tr -d '\n' | tr -d ' ')

# Run reads through a Centrifuge database
centrifuge --mm -q \
    -x /centrifuge_2023_05_16 \
    -U $data \
    --report-file $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.report \
    -S $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    --threads 12 \
    --min-hitlen 16 -k 1

# Filter out Mycobacterium reads from each paired fastq file
python3 /filterMTBreads_Kra2.py \
    $PWD/ONT_ANALYSIS/important_files/${barcode}.kraken2.txt Mycobacterium \
    $PWD/ONT_ANALYSIS/temp_files1/${barcode}.MTBReadID \
    $PWD/ONT_ANALYSIS/temp_files1/${barcode}.MTBreadNum

zcat ${files} | seqtk subseq - ${run}_${barcode}.MTBReadID > $PWD/ONT_ANALYSIS/fastq/${run}_${barcode}_Myco_reads.fastq

# Map reads to the reference and convert to SAM, BAM, and sorted BAM files
minimap2 -ax map-ont ${reference} $PWD/ONT_ANALYSIS/fastq/${run}_${barcode}_Myco_reads.fastq -t 4 > $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sam

samtools view -bq 50 $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sam --threads 4 > $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.bam

samtools sort $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.bam -o $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam --threads 4

samtools index $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam -@ 1

# Create a text file with the coverage at each position
samtools depth -aa $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam > $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_cov.txt

# Generate a Centrifuge report
centrifuge-kreport \
    -x $database \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    --min-score 150 > $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_report.txt

# Gather important files and produce output
python3 /report_resultsV1.py \
    ${run} \
    ${barcode} \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_report.txt \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_cov.txt \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    $PWD/ONT_ANALYSIS/output/${run}_${barcode}_masked.csv

# Produce a coverage plot
samtools faidx ${reference}
samtools index $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam
pysamstats -f ${reference} --type variation $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam --min-baseq 0 --pad > $PWD/ONT_ANALYSIS/output/${run}${barcode}control_samStatsUpdate.txt

# Produce a VCF file
python3 /plotBam_V1.py $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam $PWD/ONT_ANALYSIS/output/Sample_plot_${barcode}.pdf
