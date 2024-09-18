#!/usr/bin/env python3
import sys
import pandas as pd
import csv

# Define a range for coverage calculations
coverage_thresholds = range(10, 1010, 10)

def filter_mycobacterium(row):
  """Flag rows where the species name contains 'Mycobacterium'."""
  row['Myco'] = 'Mycobacterium' in row['name']
  return row

# Define headers for the output CSV
headers = [
  'Run name', 'Barcode', 'reference ID', 'total reads (n)', 'Mean read length',
  'human reads (n)', '% total reads human', 'Top Mycobacterium species name',
  'Mycobacteria reads (n) for the species identified (16s rRNA reads masked)',
  'Total reads Mycobacteria', '% total reads Mycobacteria', 'unclassified',
  '% unclassified', 'Other reads (n) not human or mycobacteria', '% other reads',
  'avDepth', 'Cov1rate', 'Cov5rate'
]

def get_species(centrifuge_report, human_filter, human_extra):
  """Process the centrifuge report to extract species-related metrics."""
  column_names = ['percentage', 'reads', 'reads_taxon', 'taxon', 'taxid', 'name']
  df = pd.read_csv(centrifuge_report, sep='\t', names=column_names, skipinitialspace=True)
  df['name'] = df['name'].str.strip()

  # Read human filter data
  df_filt = pd.read_csv(human_filter, sep='\t', names=column_names, skipinitialspace=True)
  df_filt['name'] = df_filt['name'].str.strip()

  # Initialize a dictionary for results
  results = {}

  # Total reads (n)
  total_reads = df_filt['reads'][0] + df_filt['reads'][1]
  results['total reads (n)'] = int(total_reads)

  # Unclassified reads
  unclassified = df_filt['reads'][0]
  results['unclassified'] = unclassified
  results['% unclassified'] = (unclassified / total_reads) * 100

  # Human reads (n)
  human_reads = df_filt[df_filt['taxid'] == 9606]['reads'].max() + int(human_extra)
  results['human reads (n)'] = int(human_reads)

  # % total reads human
  percent_human_reads = (float(human_reads) / float(total_reads)) * 100
  results['% total reads human'] = percent_human_reads

  # Filter for Mycobacterium species
  df = df.apply(filter_mycobacterium, axis=1)
  myco_species = df[df['Myco'] & (df['taxon'] == 'S')]
  top_myco_reads = myco_species['reads'].max()
  top_myco_name = myco_species[myco_species['reads'] == top_myco_reads]['name'].unique()[0]
  results['Mycobacteria reads (n) for the species identified (16s rRNA reads masked)'] = int(top_myco_reads)
  results['Top Mycobacterium species name'] = top_myco_name

  # % total reads Mycobacteria
  myco_reads = df[df['taxid'] == 1763]['reads'].max()
  percent_myco_reads = (myco_reads / total_reads) * 100
  results['Total reads Mycobacteria'] = float(myco_reads)
  results['% total reads Mycobacteria'] = float(percent_myco_reads)

  # Other reads (n) not human or mycobacteria
  other_reads = total_reads - (human_reads + myco_reads + unclassified)
  percent_other_reads = (other_reads / total_reads) * 100
  results['Other reads (n) not human or mycobacteria'] = int(other_reads)
  results['% other reads'] = float(percent_other_reads)

  return results

def get_coverage(coverage_file):
  """Process the coverage file to extract coverage-related metrics."""
  df = pd.read_csv(coverage_file, names=['chrom', 'pos', 'depth'], sep='\t')

  # Average depth of coverage
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

  # Initialize results dictionary
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
      if coverage_rate == 0:
          break

  coverage_results.update(additional_coverage_rates)
  return coverage_results

def get_mean_length(results_file):
  """Calculate the mean read length from the results file."""
  with open(results_file, 'r') as f:
      reader = csv.reader(f, delimiter='\t')
      fourth_column = [row[3] for row in reader]

  # Extract the numbers on either side of the pipe and calculate their sum
  numbers = [int(num1) + int(num2) for num1, num2 in (element.split('|') for element in fourth_column)]

  # Calculate the mean of the sums
  mean_length = round(sum(numbers) / len(numbers), 2)
  return mean_length

###### Inputs ########

run_name = sys.argv[1]
barcode = sys.argv[2]
centrifuge_report = sys.argv[3]
depth_coverage = sys.argv[4]
cent_results = sys.argv[5]
human_filter = sys.argv[6]
output_report = sys.argv[7]
human_extra = sys.argv[8]

# Calculate mean read length
mean_length = get_mean_length(cent_results)

###### Run functions
species_results = get_species(centrifuge_report, human_filter, human_extra)
coverage_results = get_coverage(depth_coverage)

# Combine results
combined_results = {**species_results, **coverage_results}
combined_results['Run name'] = run_name
combined_results['Barcode'] = barcode
combined_results['Mean read length'] = mean_length

# Output results to CSV
df_output = pd.DataFrame([combined_results])
df_output = df_output[headers]
df_output.to_csv(output_report, index=False)