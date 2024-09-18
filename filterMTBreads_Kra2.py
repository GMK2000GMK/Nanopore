#!/usr/bin/env python3
import sys
import pandas as pd
from ete3 import NCBITaxa

# Initialize NCBITaxa for taxonomy operations
ncbi = NCBITaxa()

# Define column names for the Kraken file
column_names = ['classifier', 'readID', 'taxID']

def get_kraken_data(kraken_file):
  """Read the Kraken file and return a DataFrame with relevant columns."""
  df = pd.read_csv(kraken_file, sep='\t', header=None, names=column_names, usecols=[0, 1, 2])
  df['taxID'] = df['taxID'].map(int)  # Ensure taxID is an integer
  return df

def get_taxids(taxon):
  """Retrieve all descendant taxIDs for a given taxon."""
  descendants = ncbi.get_descendant_taxa(taxon, collapse_subspecies=False, intermediate_nodes=True)
  name_to_taxid = ncbi.get_name_translator([taxon])
  descendants.append(int(name_to_taxid[taxon][0]))  # Add the taxon itself
  return descendants

def filter_dataframe(df, taxids):
  """Filter the DataFrame to include only rows with taxIDs in the given list."""
  filtered_df = df[df['taxID'].isin(taxids)]
  unique_reads = filtered_df[['readID']].drop_duplicates()  # Remove duplicate readIDs
  return unique_reads

##### Inputs #####

kraken_file = sys.argv[1]  # Path to the Kraken file
taxon = sys.argv[2]        # Taxon name to filter by
output_reads = sys.argv[3] # Output file for filtered readIDs
output_report = sys.argv[4] # Output file for the report

##### Run ######

# Process the Kraken file to get a DataFrame
kraken_df = get_kraken_data(kraken_file)

# Get taxIDs for the specified taxon
taxids = get_taxids(taxon)

# Filter the DataFrame to get unique readIDs for the taxon
filtered_reads = filter_dataframe(kraken_df, taxids)

# Save the filtered readIDs to a CSV file
filtered_reads.to_csv(output_reads, index=False)

# Create a report with the number of rows in the original DataFrame
report_content = 'nrow.data.frame.centReadID..\n{0}\n'.format(len(kraken_df))

# Write the report to a file
with open(output_report, 'wt') as report_file:
  report_file.write(report_content)