#!/bin/bash

# Set the names of your input files here
genome_fasta="/home/muhammedk/masters/database/complete_seq/genome.fasta"
aminoacid_fasta="/home/muhammedk/masters/database/complete_seq/amino.fasta"

# Step 0: Extract taxon IDs from genome and amino acid FASTA files
grep -E "^>" "$genome_fasta" | awk -F'|' '{print \$3}' > genome_taxon_ids.txt
grep -E "^>" "$aminoacid_fasta" | awk -F'|' '{print \$3}' > aminoacid_taxon_ids.txt

# Combine the taxon IDs from both files into one file
cat genome_taxon_ids.txt aminoacid_taxon_ids.txt > taxon_ids.txt

# Step 1: Download taxdump files
wget -q https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xf taxdump.tar.gz names.dmp nodes.dmp delnodes.dmp merged.dmp

# Step 2: Generate the conversion table
taxonkit lineage --data-dir . < taxon_ids.txt > conversion_table.txt

# Step 3: Format the conversion table
sed -i 's/;/\t/g' conversion_table.txt
awk -F '\t' '{print \$1"\t"\$3"\t"\$2}' conversion_table.txt > formatted_conversion_table.txt

# Step 4: Build Centrifuge index using the formatted conversion table
# Make sure to replace "/path/to/sequences.fasta" with the actual path to your FASTA file
centrifuge-build --conversion-table formatted_conversion_table.txt --taxonomy-tree nodes.dmp --name-table names.dmp --min-length 100 -p 4 /path/to/sequences.fasta index

# Clean up temporary files
rm -rf genome_taxon_ids.txt aminoacid_taxon_ids.txt taxon_ids.txt conversion_table.txt formatted_conversion_table.txt names.dmp nodes.dmp delnodes.dmp merged.dmp taxdump.tar.gz

echo "Custom Centrifuge database creation complete."

# Example of running centrifuge-build with specific paths
centrifuge-build --conversion-table /home/muhammedk/masters/centrifudge/test1/taxid.map --taxonomy-tree /home/muhammedk/masters/centrifudge/test1/nodes.dmp --name-table /home/muhammedk/masters/centrifudge/test1/names.dmp /home/muhammedk/masters/database/complete_seq/genome.fasta /home/muhammedk/masters/database/amino.fasta ex

# Another example of running centrifuge-build
centrifuge-build --conversion-table /home/muhammedk/masters/kraken2/kraken2_database/seqid2taxid.map --taxonomy-tree /home/muhammedk/masters/centrifudge/test1/nodes.dmp --name-table /home/muhammedk/masters/centrifudge/test1/names.dmp /home/muhammedk/masters/database/complete_seq/genome.fasta /home/muhammedk/masters/database/amino.fasta ex