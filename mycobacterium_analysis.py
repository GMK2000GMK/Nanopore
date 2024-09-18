#!/usr/bin/env python3
import sys
import math
import pandas as pd

# Define a range for coverage calculations
coverage_thresholds = range(10, 1000, 10)

def filter_mycobacterium(row):
  """Flag rows where the species name starts with 'Mycobacterium'."""
  row['IsMycobacterium'] = row['name'].startswith('Mycobacterium')
  return row

def get_species(centrifuge_report):
  """Process the centrifuge report to extract species-related metrics."""
  column_names = ['percentage', 'reads', 'reads_taxon', 'taxon', 'taxid', 'name']
  df = pd.read_csv(centrifuge_report, sep='\t', names=column_names, skipinitialspace=True)
  df['name'] = df['name'].str.strip()

  # Initialize a dictionary for results
  results = {}

  # Total reads (n)
  total_reads = df['reads_taxon'].sum()
  results['total reads (n)'] = int(total_reads)

  # Human reads (n)
  human_reads = df[df['taxid'] == 9606]['reads'].max()
  results['human reads (n)'] = int(human_reads) if not math.isnan(human_reads) else 0

  # % total reads human
  percent_human_reads = (float(human_reads) / float(total_reads)) * 100
  results['% total reads human'] = percent_human_reads

  # Mycobacteria reads (n) for the species identified (16s rRNA reads masked)
  df = df.apply(filter_mycobacterium, axis=1)
  mycobacterium_species = df[df['IsMycobacterium']]
  mycobacterium_species = mycobacterium_species[mycobacterium_species['taxon'] == 'S']
  top_myco_reads = mycobacterium_species['reads'].max()
  top_myco = mycobacterium_species[mycobacterium_species['reads'] == top_myco_reads]
  top_myco_name = top_myco['name'].unique()[0]
  results['Mycobacteria reads (n) for the species identified (16s rRNA reads masked)'] = int(top_myco_reads)
  results['Top Mycobacterium species name'] = top_myco_name

  # % total reads Mycobacteria
  mycobacterium_reads = df[df['taxid'] == 1763]['reads'].max()
  percent_myco_reads = (mycobacterium_reads / total_reads) * 100
  results['% total reads Mycobacteria'] = float(percent_myco_reads)

  # Other reads (n) not human or mycobacteria
  other_reads = total_reads - (human_reads + mycobacterium_reads)
  percent_other_reads = (other_reads / total_reads) * 100
  results['Other reads (n) not human or mycobacteria'] = int(other_reads) if not math.isnan(other_reads) else 0
  results['% other reads'] = float(percent_other_reads)

  return results

def get_coverage(coverage_file):
  """Process the coverage file to extract coverage-related metrics."""
  df = pd.read_csv(coverage_file, names=['chrom', 'pos', 'depth'], sep='\t')

  # Average depth of coverage of Mycobacteria genome for the species identified
  average_depth = df['depth'].mean()

  # Cov1rate (% genome covered at least once)
  genome_length = df['pos'].max()
  coverage_1x = df[df['depth'] >= 1]
  cov1_rate = (len(coverage_1x) / genome_length) * 100

  # Cov5rate (% genome covered at least five times)
  coverage_5x = df[df['depth'] >= 5]
  cov5_rate = (len(coverage_5x) / genome_length) * 100

  # Reference ID
  reference_id = df['chrom'].unique()[0]

  coverage_results = {
      'avDepth': average_depth,
      'Cov1rate': cov1_rate,
      'Cov5rate': cov5_rate,
      'reference ID': reference_id
  }

  # Calculate additional coverage rates for specified thresholds
  additional_coverage_rates = {}
  for threshold in coverage_thresholds:
      coverage_nx = df[df['depth'] >= threshold]
      coverage_rate = int((len(coverage_nx) / genome_length) * 100)
      additional_coverage_rates[f'Cov{threshold}rate'] = coverage_rate

  coverage_results.update(additional_coverage_rates)
  return coverage_results

def get_mean_length(results_file):
  """Calculate the mean read length from the results file."""
  df = pd.read_csv(results_file, sep='\t')
  mean_length = df['queryLength'].mean()
  return mean_length

###### Inputs ########

run_name = sys.argv[1]
barcode = sys.argv[2]
centrifuge_report = sys.argv[3]
depth_coverage = sys.argv[4]
cent_results = sys.argv[5]
output_report = sys.argv[6]

# Add additional coverage headers
headers = [
  'Run name', 'Barcode', 'Top Mycobacterium species name', 'reference ID',
  'total reads (n)', 'Mean read length', 'human reads (n)', '% total reads human',
  'Mycobacteria reads (n) for the species identified (16s rRNA reads masked)',
  '% total reads Mycobacteria', 'Other reads (n) not human or mycobacteria',
  '% other reads', 'avDepth', 'Cov1rate', 'Cov5rate'
]
for threshold in coverage_thresholds:
  headers.append(f'Cov{threshold}rate')

###### Run functions
species_results = get_species(centrifuge_report)
coverage_results = get_coverage(depth_coverage)
mean_length = get_mean_length(cent_results)

# Combine results
combined_results = {**species_results, **coverage_results}
combined_results['Run name'] = run_name
combined_results['Barcode'] = barcode
combined_results['Mean read length'] = mean_length

# Output
df_output = pd.DataFrame([combined_results])
df_output = df_output[headers]
df_output.to_csv(output_report, index=False)