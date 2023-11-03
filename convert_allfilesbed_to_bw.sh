#!/bin/bash

# Path to your chrom.sizes file
chrom_sizes="./chrom.sizes"

# Directory to store the output .bw files
output_dir="./bw_files"
mkdir -p $output_dir

# Convert each bedGraph file to BigWig
for bdg_file in ./*/*.bdg; do
  # Define the output filename based on the bdg file
  output_bw="${output_dir}/$(basename ${bdg_file%.bdg}).bw"
  
  # Convert the bedGraph to BigWig
  bedGraphToBigWig ${bdg_file} ${chrom_sizes} ${output_bw}
done

echo "Conversion complete!"
