# Mycobacterium Analysis Workflow
This repository contains a set of scripts and a workflow designed to analyze nanopore sequencing data, specifically focusing on Mycobacterium species. The workflow processes raw sequencing data, filters for Mycobacterium reads, maps them to a reference genome, and generates various reports and plots.

## Overview
The workflow is designed to:

Classify nanopore sequencing reads using a Centrifuge database.
Filter reads to retain only those belonging to Mycobacterium species.
Map filtered reads to a reference genome and generate coverage data.
Produce a variety of outputs including coverage plots, VCF files, and summary reports.
## Requirements
Python 3.7,
Bash,
pysam,
pandas,
matplotlib,
ete3,
Samtools,
Minimap2,
Centrifuge,
Seqtk

## Installation
Clone the repository:
```bash
git clone https://github.com/yourusername/mycobacterium-analysis-workflow.git
cd mycobacterium-analysis-workflow
```
## Install the required Python packages:
```bash
pip install pysam pandas matplotlib ete3
```
Ensure that Samtools, Minimap2, Centrifuge, and Seqtk are installed and available in your PATH.
Usage
## Prepare your input files
- A reference genome file.
- Fastq files containing nanopore sequencing reads.
## Run the workflow using the provided shell script:
```bash
./Nanopore_workflow.sh
```
The script will guide you through setting up the necessary directories and running the analysis.
## Scripts
plotBam_V1.py: Generates coverage plots from BAM files.

report_resultsV1.py: Processes coverage and species data to generate a summary report.

filterMTBreads_Kra2.py: Filters Mycobacterium reads from Kraken output.

Nanopore_workflow.sh: Main workflow script to automate the analysis process.

mycobacterium_analysis.py: Analyzes centrifuge reports and coverage data to extract species metrics.


## Outputs

Coverage plots: Visual representation of read coverage across the reference genome.

VCF files: Variant call format files for detected variants.

CSV reports: Summarized data including read counts, species identification, and coverage statistics.
