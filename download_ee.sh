#!/bin/bash

# Base URL
base_url="https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/entries"

# Path to the CSV file
csv_file="jiang.csv"

# Output directory for downloaded files
output_dir="jiang_frag"

# Extract the ID column, skipping the header
tail -n +2 "$csv_file" | cut -d',' -f1 | while read id; do
  # Strip whitespace
  id=$(echo "$id" | xargs)

  # Construct the download URL and output file path
  download_url="$base_url/$id/hg38/${id}.hg38.frag.tsv.bgz"
  file_path="$output_dir/${id}.hg38.frag.tsv.bgz"

  # Check if the file already exists
  if [ -f "$file_path" ]; then
    echo "$file_path already exists, skipping download."
  else
    echo "Downloading $file_path..."
    wget -O "$file_path" "$download_url"
  fi
done