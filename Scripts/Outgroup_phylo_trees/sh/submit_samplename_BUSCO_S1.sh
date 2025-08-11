#!/bin/bash

REGION=$1
BASE_DIR=$2
SAMPLE="S_868442_S1"

# Converting sample names within fasta files to NSWIDs
for file in ${BASE_DIR}/data/BUSCO/Rerun/data/fasta_BUSCO_trim/*.fa; do
  if [[ $file == *"${SAMPLE}"* ]]; then
    sed -i "s/${REGION}/${SAMPLE}/" "$file"
    echo "Renamed $file with "${SAMPLE}""
  fi
done
