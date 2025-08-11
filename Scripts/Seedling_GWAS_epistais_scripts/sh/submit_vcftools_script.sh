#!/bin/bash

REGION=$1
BASE_DIR=$2

LOG_FILE_1="${BASE_DIR}/data/log_filtvcf_DP6_GQ20.txt"
LOG_FILE_2="${BASE_DIR}/data/log_filtvcf_indels_biallelic.txt"
LOG_FILE_3="${BASE_DIR}/data/log_filtvcf_depth.txt"

# Remove indels and keep biallelic
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/rawvcf/mq.call.${REGION}.vcf.gz --keep ${BASE_DIR}/data/meta/NSWID_library_list.txt --remove-indels --min-alleles 2 --max-alleles 2  --recode --recode-INFO-all --stdout 2>> $LOG_FILE_2 | gzip -c > ${BASE_DIR}/data/filtvcf_indels_biallelic/mq.call.${REGION}_gtfilt.vcf.gz

# Calculating depth metrics
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_indels_biallelic/mq.call.${REGION}_gtfilt.vcf.gz --depth --out ${BASE_DIR}/data/qualvcf/mq.depth.${REGION}.out 2> $LOG_FILE_3