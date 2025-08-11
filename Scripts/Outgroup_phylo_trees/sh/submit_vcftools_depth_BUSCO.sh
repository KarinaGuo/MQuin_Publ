#!/bin/bash

REGION=$1
BASE_DIR=$2


# Calculating depth per indiv mean
/usr/bin/vcftools --gzvcf ${BASE_DIR}/data/filtvcf_GQ20_BUSCO/mq.call.${REGION}_gtfilt.vcf.gz --depth --out ${BASE_DIR}/data/BUSCO/evalmetric_GQ20_BUSCO/evalmetric_GQ20_${REGION}_BUSCO.out
