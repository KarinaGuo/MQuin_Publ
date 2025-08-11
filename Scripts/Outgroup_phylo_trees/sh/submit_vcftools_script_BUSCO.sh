#!/bin/bash

REGION=$1
BASE_DIR=$2


# GQ20 DP4
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/rawvcf_BUSCO/mq.call.${REGION}.vcf.gz --keep ${BASE_DIR}/data/meta/geno_out_library_list_subset.txt --remove-indels --min-alleles 2 --max-alleles 2 --minDP 4 --minGQ 20 --recode --recode-INFO-all --stdout | bgzip -c > ${BASE_DIR}/data/filtvcf_DP4_GQ20_BUSCO/mq.call.${REGION}_gtfilt.vcf.gz


