#!/bin/bash
REGION=$1
BASE_DIR=$2

LOG_FILE=${BASE_DIR}/data/log_MQ.txt

echo BASE_DIR set to ${BASE_DIR}

echo "/data/genomics/apps/bcftools-1.17/bcftools mpileup --min-MQ $MQ --annotate FORMAT/DP,FORMAT/AD -Ou -f /stage/jgb/mq1/ref/Mqui_v1_hapA.fasta -b /home/jgb/mqgwas/data/meta/bam.list.mq.all -r ${REGION} 2>> ${LOG_FILE} | /data/genomics/apps/bcftools-1.17/bcftools call -mv -f GQ -Oz -o ${BASE_DIR}/data/rawvcf/mq.call.${REGION}.vcf.gz" >> ${LOG_FILE}

/data/genomics/apps/bcftools-1.17/bcftools mpileup --min-MQ 15 --annotate FORMAT/DP,FORMAT/AD -Ou -f /stage/jgb/mq1/ref/Mqui_v1_hapA.fasta -b /home/jgb/mqgwas/data/meta/bam.list.mq.all -r ${REGION} 2>>${LOG_FILE} | /data/genomics/apps/bcftools-1.17/bcftools call -mv -f GQ -Oz -o ${BASE_DIR}/data/rawvcf/mq.call.${REGION}.vcf.gz

# add -f GQ to call? call -Oz


