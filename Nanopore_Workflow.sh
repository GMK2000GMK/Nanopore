# # TB_Workflow
# This workflow is made to analyse nanopore data in Linux, . You are expected to provide a reference sequence and fastq file. The fastq files is classified using a Centrifuge database. The results are then filtered for mycobacterium read ids. The filtered mycobacterium fastq file are mapped to your selected reference and converted into a SAM file, then BAM file and finally a sorted BAM file using Samtools. The coverage at each position is obtained. Once these processes are complete the important files are ran through a script which outputs a CSV file with the calculated results from the run. Finally from the sorted BAM file a Coverage plot and VCF file are produced.
run=''
barcode=''
ref=
db=
# Make necessary folders
# 
cd
mkdir ONT_ANALYSIS
cd ONT_ANALYSIS
mkdir fastq
mkdir temp_files1
mkdir important_files
mkdir output


# Assign the reference, barcode and run name to create a combined single fastq file
# 
echo $barcode
echo $run
files=/mnt/data/analysis/muhammedk/outbreak_bern/fastq_all/
data=`ls -m ${files} | tr -d '\n' | tr -d ' '`


# Run reads through a Centrifuge database

centrifuge --mm -q \
-x /mnt/data/analysis/muhammedk/MSC4/step2/stepB/databases/centrifuge_2023_05_16/ex \
    -U $data \
    --report-file $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.report \
    -S $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    --threads 12 \
    --min-hitlen 16 -k 1


# Filter out Mycobacterium reads from the each paired fastq file

python3 /mnt/data/analysis/muhammedk/Nanopore/tb_workflow-main-Nanopore/Nanopore/Important_scripts/filterMTBreads_Kra2.py\
    $PWD/ONT_ANALYSIS/important_files/${barcode}.kraken2.txt Mycobacterium\
    $PWD/ONT_ANALYSIS/temp_files1/${barcode}.MTBReadID\
    $PWD/ONT_ANALYSIS/temp_files1/${barcode}.MTBreadNum

zcat ${files} | seqtk subseq - ${run}_${barcode}.MTBReadID > $PWD/ONT_ANALYSIS/fastq/${run}_${barcode}_Myco_reads.fastq


# Assign read variables to filtered paired fastq files

minimap2 -ax map-ont ${ref} $PWD/ONT_ANALYSIS/fastq/${run}_${barcode}_Myco_reads.fastq -t 4 > $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sam

# convert sam file to bam file and filter our any read wih mapq lower than 50
samtools view -bq 50 $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sam --threads 4 > $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.bam

# sort the bam file with 4 cores
samtools sort $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.bam -o $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam --threads 4

# index the sorted bam file file
samtools index $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam -@ 1

# create a text file witg the coverage at each position
samtools depth -aa $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam > $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_cov.txt



centrifuge-kreport \
    -x $db \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    --min-score 150 > $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_report.txt

    
# Gather important files and produce output

python3 /mnt/data/analysis/muhammedk/Nanopore/tb_workflow-main-Nanopore/Nanopore/Important_scripts/report_resultsV1.py \
    ${run} \
    ${barcode} \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_report.txt \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_cov.txt \
    $PWD/ONT_ANALYSIS/important_files/${run}_${barcode}_Cent_K1.results \
    $PWD/ONT_ANALYSIS/output/${run}_${barcode}_masked.csv

# Produce a coverage plot

samtools faidx ${refs}
samtools index $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam
pysamstats -f ${refs} --type variation $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam --min-baseq 0 --pad > $PWD/ONT_ANALYSIS/output/${run}${barcode}control_samStatsUpdate.txt


# Produce a VCF file

python3 /mnt/data/analysis/muhammedk/Nanopore/tb_workflow-main-Nanopore/Nanopore/Important_scripts/plotBam_V1.py $PWD/ONT_ANALYSIS/temp_files1/${run}_${barcode}.sorted.bam $PWD/ONT_ANALYSIS/output\Sample_plot_${barcode}.pdf
