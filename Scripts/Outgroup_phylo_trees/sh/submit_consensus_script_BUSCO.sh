#!/bin/bash

REGION=$1
BASE_DIR=$2


# Converting VCF to BCF to allow consensus

/data/genomics/apps/bcftools-1.17/bcftools view -O b -o ${BASE_DIR}/data/filtbcf_DP4_GQ20_BUSCO/mq.call.${REGION}_gtfilt.bcf ${BASE_DIR}/data/filtvcf_DP4_GQ20_BUSCO/mq.call.${REGION}_gtfilt.vcf.gz

# Index

/data/genomics/apps/bcftools-1.17/bcftools index ${BASE_DIR}/data/filtbcf_DP4_GQ20_BUSCO/mq.call.${REGION}_gtfilt.bcf

# Convert to fasta by region and individual; trimming generated fasta file to region of interest

while read SAMPLE; do 
  cat /stage/jgb/mq1/ref/Mqui_v1_hapA.fasta | /data/genomics/apps/bcftools-1.17/bcftools consensus -s $SAMPLE --missing "N" --mark-del "N" --mark-ins lc -I ${BASE_DIR}/data/filtbcf_DP4_GQ20_BUSCO/mq.call.${REGION}_gtfilt.bcf > ${BASE_DIR}/data/fasta_BUSCO/mq.call.${REGION}_${SAMPLE}.fa
  echo "Consensus generated for ${REGION} ; ${SAMPLE}"
  /usr/bin/samtools faidx ${BASE_DIR}/data/fasta_BUSCO/mq.call.${REGION}_${SAMPLE}.fa ${REGION} -o ${BASE_DIR}/data/fasta_BUSCO_trim/mq.trim.${REGION}_${SAMPLE}.fa
  rm ${BASE_DIR}/data/fasta_BUSCO/mq.call.${REGION}_${SAMPLE}.fa
  echo "Consensus trimmed for ${REGION} ; ${SAMPLE}"
done < ${BASE_DIR}/data/meta/geno_out_library_list_subset.txt