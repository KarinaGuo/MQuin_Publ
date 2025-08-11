#!/bin/bash

REGION=$1
BASE_DIR=$2

# Converting sample names within fasta files to NSWIDs

while read SAMPLE; do
  for file in ${BASE_DIR}/data/fasta_BUSCO_trim/*.fa; do
    if [[ $file == *"${SAMPLE}"* ]]; then
      sed -i "1s/^>.*/>${SAMPLE}/" "$file"
      echo "Renamed $file with "${SAMPLE}""
    fi
  done
done < ${BASE_DIR}/data/meta/geno_out_library_list_subset.txt
