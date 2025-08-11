# BLASTing
BASE_DIR=/home/karina/mqgwas/snpeff/g93tog98/
cd /home/karina/mqgwas/snpeff/g93tog98/

while IFS=$'\t' read -r gene REGION; do
 echo ${gene} ${REGION}
 /usr/bin/samtools faidx /stage/jgb/mq1/ref/Mqui_v1_hapA.fasta ${REGION} -o fasta/mq.trim.${gene}_reference_HapA.fa

 for hap in hapA hapB; do
  /data/genomics/apps/blast/ncbi-blast-2.14.0+/bin/blastn -db /home/jgb/mqgwas/hits/dbs/Mqui_v1_${hap}.fasta -query fasta/mq.trim.${gene}_reference_HapA.fa -out blast/mq.blast.${gene}_reference_${hap}.fa -evalue 1e-6 -outfmt 6
 done

done < gene_region

# concatenate to one file

while IFS=$'\t' read -r gene REGION; do
 touch blast/concat/mq.blast.concat.${gene}.txt
 for hap in hapA hapB; do
  cat blast/mq.blast.${gene}_reference_${hap}.fa >> blast/concat/mq.blast.concat.${gene}.txt
 done
done < gene_region

while IFS=$'\t' read -r gene REGION; do
 awk -F'\t' '$2 == "MqB_CHR08" && $12 > 3000' blast/concat/mq.blast.concat.${gene}.txt > blast/concat/mq.top_HapBres_blast.${gene}.txt
done < gene_region

## MAF via anchorwave
## for the gene region find the corresponding hapB region
grep -E '^[as]' /home/karina/mqgwas/snpeff/g93tog98/Mqui_hapB.f.slice.maf.txt | cut -c1-100 > headers_maf.txt ## Extract any lines that start with "a" or "s" & the first 100 characters of the line
