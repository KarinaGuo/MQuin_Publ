#!/bin/bash

REGION=$1
BASE_DIR=$2

LOG_FILE_4="${BASE_DIR}/data/log_exlowdepth.txt"
LOG_FILE_5="${BASE_DIR}/data/log_exmissing.txt"
LOG_FILE_1="${BASE_DIR}/data/log_filtvcf_DP6_GQ20.txt"
LOG_FILE_2="${BASE_DIR}/data/log_filtvcf_indels_biallelic.txt"
LOG_FILE_3="${BASE_DIR}/data/log_filtvcf_depth.txt"


# Remove individuals with a weighted mean coverage depth less than 3
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_indels_biallelic/mq.call.${REGION}_gtfilt.vcf.gz --keep ${BASE_DIR}/data/meta/samples_gt_depth_thresh_indel-biallelic.txt --recode --recode-INFO-all --stdout 2> $LOG_FILE_4 | gzip -c > ${BASE_DIR}/data/filtvcf_exlowdepth/mq.call.${REGION}_gtdepth.vcf.gz 


# Mask sites with lt DP>6 GQ 20
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_exlowdepth/mq.call.${REGION}_gtdepth.vcf.gz  --minDP 6 --minGQ 20 --recode --recode-INFO-all --stdout 2> $LOG_FILE_1 | gzip -c > ${BASE_DIR}/data/filtvcf_DP6_GQ20/mq.call.${REGION}_gtmask.vcf.gz 


# Filter on missingness
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_DP6_GQ20/mq.call.${REGION}_gtmask.vcf.gz --max-missing-count 200 --min-alleles 2 --recode --recode-INFO-all --stdout 2> $LOG_FILE_5 | gzip -c > ${BASE_DIR}/data/filtvcf_missing/mq.call.${REGION}_gtmask.vcf.gz

