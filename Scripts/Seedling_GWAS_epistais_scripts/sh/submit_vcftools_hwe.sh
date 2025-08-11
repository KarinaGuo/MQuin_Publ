#!/bin/bash

REGION=$1
BASE_DIR=$2
LOG_FILE=${BASE_DIR}/data/log_hwe.txt

# Filter on HWE
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_missing/mq.call.${REGION}_gtmask.vcf.gz --hwe 1e-4 --recode --recode-INFO-all --stdout 2>> $LOG_FILE | gzip -c > ${BASE_DIR}/data/filtvcf_hwe/mq.call.${REGION}_hwe.vcf.gz