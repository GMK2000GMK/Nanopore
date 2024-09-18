#!/usr/bin/env python3
import pysam
import sys
import matplotlib as mpl
mpl.use('Agg')  # Use non-interactive backend for matplotlib
import matplotlib.pyplot as plt

def get_stats(bam_file):
  """Extract basic statistics from a BAM file using pysam."""
  keys = ['chrom', 'length', 'mapped', 'unmapped']
  stats = dict(zip(keys, pysam.idxstats(bam_file).split('\n')[0].split('\t')))
  return stats

def get_coverage(bam_file):
  """Calculate coverage depth across a reference sequence in a BAM file."""
  bam_stats = get_stats(bam_file)
  coverage = []
  samfile = pysam.AlignmentFile(bam_file, "rb")
  current_position = 0

  # Iterate over each position in the reference sequence
  for pileupcolumn in samfile.pileup(bam_stats['chrom'], 0, int(bam_stats['length'])):
      # Fill in zero coverage for positions with no reads
      if pileupcolumn.pos > current_position:
          coverage.extend([0] * (pileupcolumn.pos - current_position))
          current_position = pileupcolumn.pos

      # Append coverage depth for the current position
      coverage.append(pileupcolumn.n)
      current_position += 1

  return coverage, bam_stats

def smooth(data, window_size):
  """Smooth the coverage data using a sliding window average."""
  data_length = len(data)
  smoothed_data = [sum(data[i:i+window_size]) / float(window_size) for i in range(0, data_length - window_size, int(data_length / 100.0))]
  return smoothed_data

def plot_coverage(coverage, bam_stats):
  """Plot the coverage depth across the reference sequence."""
  reference_name, sample_id = extract_name(sys.argv[1])
  # Optionally smooth the coverage data
  # smoothed_coverage = smooth(coverage, 100)
  smoothed_coverage = coverage

  # Set plot limits and labels
  ymax = max(smoothed_coverage)
  plt.ylim(0, ymax * 1.05)
  plt.plot(smoothed_coverage)
  plt.title("Sample {0}, ref {1}".format(sample_id, bam_stats['chrom']))
  plt.xlabel("Position on reference")
  plt.ylabel("Coverage depth (bp)")
  plt.margins(y=0)
  plt.savefig(sys.argv[2])  # Save the plot to a file

def extract_name(file_path):
  """Extract the reference name and sample ID from the file path."""
  path_parts = file_path.split('/')
  reference_name = path_parts[-1].replace(".sorted.bam", "")
  sample_id = path_parts[-1].split('.')[0]
  return reference_name, sample_id

# Main execution
coverage_data, bam_stats = get_coverage(sys.argv[1])
print('Ref length', 'total bases', 'average depth')
print(len(coverage_data), sum(coverage_data), sum(coverage_data) / float(len(coverage_data)))
plot_coverage(coverage_data, bam_stats)