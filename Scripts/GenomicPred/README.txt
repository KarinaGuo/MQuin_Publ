Protocol:
* Extract sigSNP locs from GWAS folder
* Run vcf to keep positions and extract genotype
* Convert the output vcf to tassel

## Extracting Chrom and Pos
/data/genomics/apps/bcftools-1.17/bcftools query -f '%CHROM %POS\n' "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf" > /home/karina/GAPIT_29_POS.txt
/data/genomics/apps/bcftools-1.17/bcftools query -f '%CHROM %POS\n' "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf" > /home/karina/GAPIT_16_POS.txt

/data/genomics/apps/bcftools-1.17/bcftools query -f '%CHROM %POS\n' "/home/karina/mqgwas/iter_6/data/iter_4_2_thin5000_GAPIT_DArT_BGSNPs.recode.vcf" > /home/karina/GAPIT_BGSNPs_POS.txt

######################################################################################
GAPIT_1 - Used P<0.01 from Iter_2 gwas_5 

/usr/bin/vcftools --vcf "/recer1/karina/seedlings/iter_2/gwas_5/Mq_filt_cat.beagle.vcf" --positions "/home/karina/mqgwas/iter_5/data/meta/CHR08_sig_SNPS_locs_01.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_5/data/filtvcf_sigsnps/MQ_sig-sites_01_beagle.vcf.gz # Didn't use iter_2/gwas_6's beagle file as that had been filtered for hwe, removing some sites that were present in the unimputed dataset

#Genotypes for sigsites were then retrieved for the imputed dataset

/usr/bin/vcftools --gzvcf /home/karina/mqgwas/iter_5/data/filtvcf_sigsnps/MQ_sig-sites_01_beagle.vcf.gz --extract-FORMAT-info GT --out /home/karina/mqgwas/iter_5/data/qualvcf/evalmetric_Mq_GT_sigSNP_01_beagle.out


# Convert vcf to hapmap
/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf /home/karina/mqgwas/iter_6/data/MQ_sig-sites_01.vcf -export /home/karina/mqgwas/iter_6/data/MQ_sig-sites_01_hapmap -exportType Hapmap

######################################################################################
# GAPIT_2 - Used Iter_4 gwas_2 Top summary 85 sigSNPs (all normal filters + MQ=15)
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_2_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_2_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_2_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_2_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_2_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_2_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_2_sigSNP.hapmap -exportType Hapmap

######################################################################################
# GAPIT_3 - Used Iter_4 gwas Top P-val<0.0001 211 sigSNPs (all normal filters + MQ=15 + HWE1-e4)
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filthwe_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_3_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_3_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_3_sigSNP.vcf.gz 

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_3_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_3_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_3_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_3_sigSNP.hapmap -exportType Hapmap

######################################################################################
# GAPIT_4 - Used Iter_4 gwas_2 Top P-val<0.0001 238 sigSNPs (all normal filters + MQ=15) (note, only 7 SNPs lesser than the gwas run with no MQ=15 filter
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_4_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_4_sigSNP.hapmap -exportType Hapmap

######################################################################################
# GAPIT_5 - Used Iter_4 gwas_2 Top 1000 sigSNPs (all normal filters + MQ=15) 
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_5_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_5_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_5_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_5_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_5_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_5_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_5_sigSNP.hapmap -exportType Hapmap

######################################################################################
# GAPIT_6 - Used Iter 2 gwas_2 summary 95 SNPs (all normal filters) and removed problematic SNPs --> 70 SNPs
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_6_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_6_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_6_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_6_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_6_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_6_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_6_sigSNP.hapmap -exportType Hapmap

######################################################################################
# GAPIT_7 - Used Iter 2 gwas_2 summary 95 SNPs (all normal filters) then removed the same 25 problematic SNPs based on depth, HWE, blast res, MAF --> 290 SNPs
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_7_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_7_sigSNP.hapmap -exportType Hapmap


###################################################################################
# GAPIT_8
# Thinning the dataset - 700 SNPs randomly distributed across the dataset (1 SNP in 1000BP window) -> using that to generate a kinship matrix -> using kinship matrix in GAPIT

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 20000 --recode --recode-INFO-all --stdout | gzip -c > /home/karina/mqgwas/iter_6/data/GAPIT_8/Mq_vcf_thinnedSNPs.vcf.gz
# After filtering, kept 11294 out of a possible 33269031 Sites

gzip -c /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_4_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_8/GAPIT_4_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz Mq_vcf_thinnedSNPs.vcf.gz -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_8/GAPIT_8_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_8/GAPIT_8_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_8/GAPIT_8_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_8_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_8_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_8_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_8_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_8_sigSNP.hapmap -exportType Hapmap


#####################################################################################
# GAPIT_9
# Redoing GAPIT_7 but with a HWE 1e-4 filter on SNPs
/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_7_sigSNP.vcf" --hwe 1e-4 --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_9_sigSNP

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_9_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_9_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_9_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_9_sigSNP.hapmap -exportType Hapmap


#####################################################################################
# GAPIT_10
# Redoing GAPIT_8 but with concat -a

cp *.vcf.gz /home/karina/mqgwas/iter_6/data/GAPIT_10/

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4_sigSNP.vcf.gz"
gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_10/Mq_vcf_thinnedSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4_sigSNP.vcf"> "/home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4_sigSNP.vcf.gz"
bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_10/Mq_vcf_thinnedSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_10/Mq_vcf_thinnedSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4_sigSNP.vcf.gz"
tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_10/Mq_vcf_thinnedSNPs.vcf.gz"

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_4_sigSNP.vcf.gz -D -a -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_10_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_10_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_10_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_10_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_10_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_10_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_10_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_10_sigSNP.hapmap -exportType Hapmap

#####################################################################################
# GAPIT_11
# CHR08 19 million -> 19.5 million with a MAF of 0.05 and use all of it for genomic prediction

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.05 --chr MqA_CHR08 --from-bp 19000000 --to-bp 19500000 --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP
#After filtering, kept 6494 out of a possible 33269031 Sites

rm /home/karina/mqgwas/iter_6/data/sigSNP_vcf/*log
mv /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP* /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_11_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_11_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_12
# Thinning the dataset - 1 SNP in 3000BP window -> using that in GAPIT

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 30000 --recode --recode-INFO-all --stdout | gzip -c > /home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf.gz
# After filtering, kept 7733 out of a possible 33269031 Sites

cp /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4* /home/karina/mqgwas/iter_6/data/GAPIT_12/

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf.gz"

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_4_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_12/GAPIT_12_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_12/GAPIT_12_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_12/GAPIT_12_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_12_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_12_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_12_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_12_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_12_sigSNP.hapmap -exportType Hapmap

#################################################################################
# GAPIT_13
# Thinning the dataset - 1 SNP in 6000BP window -> using that in GAPIT
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 60000 --recode --recode-INFO-all --stdout | gzip -c > /home/karina/mqgwas/iter_6/data/GAPIT_13/Mq_vcf_thinnedSNPs.vcf.gz

cp /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4* /home/karina/mqgwas/iter_6/data/GAPIT_13/

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_13/Mq_vcf_thinnedSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_13/Mq_vcf_thinnedSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_13/Mq_vcf_thinnedSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_13/Mq_vcf_thinnedSNPs.vcf.gz"

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_4_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_13/GAPIT_13_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_13/GAPIT_13_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_13/GAPIT_13_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_13_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_13_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_13_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_13_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_13_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_14
# Blasted an old Dartseq onto reference genome, extracted positions from imputed beagle vcf. 

## Taking blasted positions
docker run --rm -it -v /home/karina/mqdart/:/mqdart karinaguo/pkg_test:RRwf1and3_v19 R

blast_dart <- read.table("mqdart/haps2MqA.blast")
colnames(blast_dart) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

dart_positions <- blast_dart [blast_dart$bitscore>85,] # Filtered blast results by 0.85!

#Filter out rows of duplicate qseqid, picking only the occurrence with the highest bitscore
dart_positions$bitscore <- as.numeric(dart_positions$bitscore)
sorted_dart_positions <- dart_positions[order(-dart_positions$bitscore), ]

filtered_blast_dart <- sorted_dart_positions[!duplicated(sorted_dart_positions$qseqid), ]

filtered_blast_dart_sub <- filtered_blast_dart[,c(2,9,10)]
filtered_blast_dart_sub$sstart_sort <- ifelse(filtered_blast_dart_sub$sstart > filtered_blast_dart_sub$send, filtered_blast_dart_sub$send, filtered_blast_dart_sub$sstart)
filtered_blast_dart_sub$send_sort <- ifelse(filtered_blast_dart_sub$sstart > filtered_blast_dart_sub$send, filtered_blast_dart_sub$sstart, filtered_blast_dart_sub$send)


filtered_blast_dart_sub$positions <- paste0(filtered_blast_dart_sub$sseqid, " ", filtered_blast_dart_sub$sstart_sort, " ",  filtered_blast_dart_sub$send_sort)
filtered_blast_dart_pos <- unique(filtered_blast_dart_sub[, 6])

write.table(filtered_blast_dart_pos, file = "mqdart/dart_positions_ref.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

## 
cp "/home/karina/mqdart/dart_positions_ref.txt" /home/karina/mqgwas/iter_6/data/meta/GAPIT_14_sigSNPlocs.txt

sed 's/ /\t/g' /home/karina/mqgwas/iter_6/data/meta/GAPIT_14_sigSNPlocs.txt > /home/karina/mqgwas/iter_6/data/meta/GAPIT_14_sigSNPlocs_tmp.txt

/data/genomics/apps/bedtools/bedtools2/bin/bedtools intersect -header -wa -a "/home/karina/mqgwas/iter_4/gwas_2/Mq_filt_cat.beagle.vcf" -b "/home/karina/mqgwas/iter_6/data/meta/GAPIT_14_sigSNPlocs_tmp.txt" | gzip -c >  /home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf.gz

# /usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/gwas_2/Mq_filt_cat.beagle.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_14_sigSNPlocs.txt" --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf.gz # subsetted with bedtools

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf.gz"

cp /home/karina/mqgwas/iter_6/data/GAPIT_12/GAPIT_4_sigSNP.vcf* /home/karina/mqgwas/iter_6/data/GAPIT_14/

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_dartSNPs.vcf.gz GAPIT_4_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_14/GAPIT_14_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_14/GAPIT_14_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_14/GAPIT_14_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_14_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_14_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_14_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_14_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_14_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_15
# Grabbed the P<0.0001 SNPs of COI from a parent gwas (GWAS_5, ran with breeding value of COI and RGR). Used previous thinned dataset for kinship, note not filtered for HWE 1e-5 but MQ=15 along with others

cp /home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf* /home/karina/mqgwas/iter_6/data/GAPIT_15

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_15_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_sigSNP.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_sigSNP.vcf.gz
bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_sigSNP.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_sigSNP.vcf.gz

bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_15/Mq_vcf_thinnedSNPs.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_15/Mq_vcf_thinnedSNPs.vcf
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_15/Mq_vcf_thinnedSNPs.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_15_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_15/GAPIT_15_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_15_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_15_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_15_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_15_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_15_sigSNP.hapmap -exportType Hapmap

#################################################################################
# GAPIT_16
# Thinning the dataset - 1 SNP in 15000BP window -> using that in GAPIT (14916 SNPs)
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 15000 --recode --recode-INFO-all --stdout | gzip -c > /home/karina/mqgwas/iter_6/data/GAPIT_16/Mq_vcf_thinnedSNPs.vcf.gz # 14916 SNPs

cp /home/karina/mqgwas/iter_6/data/GAPIT_10/GAPIT_4* /home/karina/mqgwas/iter_6/data/GAPIT_16/

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_16/Mq_vcf_thinnedSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_16/Mq_vcf_thinnedSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_16/Mq_vcf_thinnedSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_16/Mq_vcf_thinnedSNPs.vcf.gz"

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_4_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_16/GAPIT_16_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_16/GAPIT_16_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_16/GAPIT_16_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_16_sigSNP.hapmap -exportType Hapmap


#################################################################################
# GAPIT_17
# Top P<0.001 SNPs from iter_4/gwas_2

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_17_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_17_sigSNP.hapmap -exportType Hapmap

#################################################################################
# GAPIT_18
# Rerun DArT gBLUP with --thin 100

/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_6/data/GAPIT_14/Mq_vcf_dartSNPs.vcf" --thin 100 --recode --recode-INFO-all --out "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf"

mv "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf.recode.vcf" "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf" 

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_18/Mq_vcf_dartSNPs_thinned.vcf.gz"

cp /home/karina/mqgwas/iter_6/data/GAPIT_12/GAPIT_4_sigSNP.vcf* /home/karina/mqgwas/iter_6/data/GAPIT_18/

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_dartSNPs_thinned.vcf.gz GAPIT_4_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_18/GAPIT_18_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_18/GAPIT_18_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_18/GAPIT_18_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_18_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_18_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_18_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_18_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_18_sigSNP.hapmap -exportType Hapmap



#################################################################################
# GAPIT_19
# Top 2000 marker SNPs of parents from breeding value gwas_5

/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_4/gwas/Mq_hwe_cat.beagle.vcf" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_19_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_19_sigSNP.vcf.gz

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_19_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_19_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_19_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_19_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_19_sigSNP.hapmap -exportType Hapmap


#################################################################################
# GAPIT_20
# Run HWE 1e-5 GWAS of parents BV -> sigSNPs -> seedling positions extracted with another HWE 1e-5 filter - Top 2000 marker SNPs of parents

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_20_sigSNPlocs.txt" --hwe 1e-5 --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_20_sigSNP.vcf.gz
# After filtering, kept 524 out of a possible 33269031 Sites, not ever sigSNP of parents were kept

# Convert vcf to hapmap
gzip -d /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_20_sigSNP.vcf.gz

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_20_sigSNP.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_20_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_20_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_20_sigSNP.hapmap -exportType Hapmap

#################################################################################
# GAPIT_21
# Run HWE 1e-5 GWAS of parents BV -> sigSNPs -> seedling positions extracted with another HWE 1e-5 filter - Top P<0.001 SNPs from iter_4/gwas_2

cp /home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf* /home/karina/mqgwas/iter_6/data/GAPIT_21

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_21_sigSNPlocs.txt" --hwe 1e-5 --recode --recode-INFO-all --stdout | gzip -c >  ~/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_sigSNP.vcf.gz
# After filtering, kept 1933 out of a possible 33269031 Sites

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_sigSNP.vcf.gz
bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_sigSNP.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat Mq_vcf_thinnedSNPs.vcf.gz GAPIT_21_sigSNP.vcf.gz -a -D -Oz -o  /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_21/GAPIT_21_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_21_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_21_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_21_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_21_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_21_sigSNP.hapmap -exportType Hapmap


#################################################################################
# GAPIT_22
# Increasing testing vs training gBLUP
# Done in R

#################################################################################
# GAPIT_23
# Generate genomic predictions for COI for the Qld samples, and the small batch of broad leaved melaleucas

mkdir /home/karina/mqgwas/iter_6/data/GAPIT_23
cp /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_17_sigSNP* /home/karina/mqgwas/iter_6/data/GAPIT_23

bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_17_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_17_sigSNP.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_17_sigSNP.vcf.gz

/usr/bin/vcftools --gzvcf "/recer1/karina/outgroups_iter_1/data/catvcf/Mq_DP4_GQ20_indels_onlybiall_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_17_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_sigSNP.vcf.gz
# After filtering, kept 2377 out of a possible 52569109 Sites

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_sigSNP.vcf.gz 
bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_sigSNP.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools merge GAPIT_23_sigSNP.vcf.gz GAPIT_17_sigSNP.vcf.gz -Oz -o /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_23_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_23_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_23_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_23_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_23_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_23_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_24
# Using GAPIT_12 thinned dataset (1 SNP in 3000BP window) and no sigSNPs

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile "/home/karina/mqgwas/iter_6/data/GAPIT_12/Mq_vcf_thinnedSNPs.vcf" -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_24_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_24_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_24_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_25 Repeating GAPIT_17 to see reliability for outlier populations

###################################################################################
# GAPIT_26
# Using GAPIT_17 sigSNPs (seedlings) on parents (HWE1e-4 MQ15) to predict BV. Cross-validating these FID BV with non-genotyped seedling BV

mkdir GAPIT_26
cd GAPIT_26

cp /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_17* .

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/parent_iter_2/data/catvcf/Mq_filtmiss_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/meta/GAPIT_17_sigSNPlocs.txt" --recode --recode-INFO-all --stdout | gzip -c >  /home/karina/mqgwas/iter_6/data/GAPIT_23/GAPIT_26_sigSNP.vcf.gz
gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_sigSNP.vcf.gz
bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_sigSNP.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_sigSNP.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools merge GAPIT_26_sigSNP.vcf.gz GAPIT_17_sigSNP.vcf.gz -Oz -o /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_26/GAPIT_26_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_26_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_26_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_26_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_26_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_26_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_27
# Using GAPIT_17 sigSNPs (seedlings) and thinned by 100

mkdir GAPIT_27
cd GAPIT_27

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf GAPIT_17_sigSNP.vcf --thin 100 --recode --recode-INFO-all --out GAPIT_17_sigSNP_thin100 #After filtering, kept 1377 out of a possible 2377 Sites

mv /home/karina/mqgwas/iter_6/data/GAPIT_27/GAPIT_17_sigSNP_thin100.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_27_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_27_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_27_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_27_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_27_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_28
# Using GAPIT_17 sigSNPs (seedlings) and thinned by 500

mkdir GAPIT_28
cd GAPIT_28

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf GAPIT_17_sigSNP.vcf --thin 500 --recode --recode-INFO-all --out GAPIT_17_sigSNP_thin500 #After filtering, kept 869 out of a possible 2377 Sites

mv /home/karina/mqgwas/iter_6/data/GAPIT_28/GAPIT_17_sigSNP_thin500.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_28_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_28_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_28_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_28_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_28_sigSNP.hapmap -exportType Hapmap

###################################################################################
# GAPIT_29
# Using GAPIT_17 sigSNPs (seedlings) and thinned by 1000

mkdir GAPIT_29
cd GAPIT_29

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf GAPIT_17_sigSNP.vcf --thin 1000 --recode --recode-INFO-all --out GAPIT_17_sigSNP_thin1000 #After filtering, kept 689 out of a possible 2377 Sites

mv /home/karina/mqgwas/iter_6/data/GAPIT_29/GAPIT_17_sigSNP_thin1000.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_29_sigSNP.hapmap -exportType Hapmap


###################################################################################
# GAPIT_30
# Using GAPIT_17 sigSNPs (seedlings) and thinned by 2000

mkdir GAPIT_30
cd GAPIT_30

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf GAPIT_17_sigSNP.vcf --thin 2000 --recode --recode-INFO-all --out GAPIT_17_sigSNP_thin2000 #After filtering, kept 689 out of a possible 2377 Sites

mv /home/karina/mqgwas/iter_6/data/GAPIT_30/GAPIT_17_sigSNP_thin2000.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_30_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_30_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_30_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_30_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_30_sigSNP.hapmap -exportType Hapmap

##################################################################################
# GAPIT_31
# Using GAPIT_17 sigSNPs (seedlings) and HWE1e-4

mkdir GAPIT_31
cd GAPIT_31

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf GAPIT_17_sigSNP.vcf --hwe 1e-4 --recode --recode-INFO-all --out GAPIT_17_sigSNP_hwe1e4 #After filtering, kept 1890 out of a possible 2377 Sites

mv /home/karina/mqgwas/iter_6/data/GAPIT_31/GAPIT_17_sigSNP_hwe1e4.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_31_sigSNP.hapmap -exportType Hapmap


##################################################################################
# GAPIT_32
# Using GAPIT_17 sigSNPs (seedlings) and HWE1e-4 thin 1000

mkdir GAPIT_32
cd GAPIT_32

cp "/recer1/karina/seedlings/iter_6/GAPIT_23/GAPIT_17_sigSNP.vcf" .

/usr/bin/vcftools --vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP.vcf --thin 200 --recode --recode-INFO-all --out GAPIT_17_sigSNP_hwe1e4_thin200 #After filtering, kept 1003 out of a possible 2377 Sites 

mv /home/karina/mqgwas/iter_6/data/GAPIT_32/GAPIT_17_sigSNP_hwe1e4_thin200.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_32_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_32_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_32_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_32_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_32_sigSNP.hapmap -exportType Hapmap

##################################################################################
# GAPIT_33
# Final GAPIT_v1? Merging iteration 16 (background SNPs) & 29 (sigSNPs)

mkdir GAPIT_33
cd GAPIT_33

cp "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf" .
bgzip -c "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_16_sigSNP.vcf.gz"
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_16_sigSNP.vcf.gz

cp "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf" .
bgzip -c "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_29_sigSNP.vcf.gz"
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_29_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat GAPIT_16_sigSNP.vcf.gz GAPIT_29_sigSNP.vcf.gz -a -D -Oz -o /home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_33_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_33_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_33/GAPIT_33_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_33_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_33_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_33_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_33_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_33_sigSNP.hapmap -exportType Hapmap


##################################################################################
# GAPIT_34
# Final GAPIT_v2? Merging iteration 16 (background SNPs) & 31 (sigSNPs)

mkdir GAPIT_34
cd GAPIT_34

cp "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf" .
bgzip -c "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_16_sigSNP.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_16_sigSNP.vcf.gz"
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_16_sigSNP.vcf.gz

cp "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP.vcf" .
bgzip -c "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_31_sigSNP.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_31_sigSNP.vcf.gz"
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_31_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat GAPIT_16_sigSNP.vcf.gz GAPIT_31_sigSNP.vcf.gz -a -D -Oz -o /home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_34_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_34_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_34/GAPIT_34_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_34_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_34_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_34_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_34_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_34_sigSNP.hapmap -exportType Hapmap

##################################################################################
# GAPIT_35
# Changing BG SNPs to 5000, thin by a 50000BP window - merging with GAPIT 31 and seeing results + checking kinship validity

# Thinning the dataset - 1 SNP in 50000BP window -> using that in GAPIT
cd /home/karina/mqgwas/iter_6/data/GAPIT_35/ 

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 50000BP --recode --recode-INFO-all --stdout | gzip -c > /home/karina/mqgwas/iter_6/data/GAPIT_35/Mq_vcf_thinnedSNPs.vcf.gz

gzip -d "/home/karina/mqgwas/iter_6/data/GAPIT_35/Mq_vcf_thinnedSNPs.vcf.gz"

bgzip -c "/home/karina/mqgwas/iter_6/data/GAPIT_35/Mq_vcf_thinnedSNPs.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_35/Mq_vcf_thinnedSNPs.vcf.gz"

tabix -p vcf "/home/karina/mqgwas/iter_6/data/GAPIT_35/Mq_vcf_thinnedSNPs.vcf.gz"

cp "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf" /home/karina/mqgwas/iter_6/data/GAPIT_35/
bgzip -c "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_29_sigSNP.vcf" > "/home/karina/mqgwas/iter_6/data/GAPIT_35/GAPIT_29_sigSNP.vcf.gz"
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_35/GAPIT_29_sigSNP.vcf.gz

/data/genomics/apps/bcftools-1.17/bcftools concat GAPIT_29_sigSNP.vcf.gz Mq_vcf_thinnedSNPs.vcf.gz -a -D -Oz -o /home/karina/mqgwas/iter_6/data/GAPIT_35/GAPIT_35_concat.vcf.gz

gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_35/GAPIT_35_concat.vcf.gz

mv /home/karina/mqgwas/iter_6/data/GAPIT_35/GAPIT_35_concat.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_35_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_35_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_35_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_35_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_35_sigSNP.hapmap -exportType Hapmap

##################################################################################
# GAPIT_36
# Rerunning GAPIT_29 but replacing the top 20 SNPs with SNPs unique to the relevant sequence

mkdir GAPIT_36
cd GAPIT_36

cp "/home/karina/mqgwas/iter_4/gwas_2/Output/Linear_Mixed_Model/COI/COI_20231025_215412/best_p-values/p_wald_COI_mod_sub_GWAS_phenotype_height_COItransformedfoursqrt.part1_Mq_filt_cat.beagle_top0.001.csv" . # Manually replace top 20 SNPs

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "/home/karina/mqgwas/iter_6/data/GAPIT_36/p_wald_COI_mod_sub_GWAS_phenotype_height_COItransformedfoursqrt.part1_Mq_filt_cat.beagle_top0.001.txt" --recode --recode-INFO-all --out GAPIT_36_unthin_sigSNP

/usr/bin/vcftools --vcf GAPIT_36_unthin_sigSNP.recode.vcf --thin 1000 --recode --recode-INFO-all --out GAPIT_36_thin1000_sigSNP #After filtering, kept 689 out of a possible 2366 Sites

bgzip -c GAPIT_36_thin1000_sigSNP.recode.vcf > GAPIT_36_thin1000_sigSNP.recode.vcf.gz
tabix -p vcf GAPIT_36_thin1000_sigSNP.recode.vcf.gz
gzip -d GAPIT_36_thin1000_sigSNP.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_36/GAPIT_36_thin1000_sigSNP.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_36_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_36_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_36_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_36_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_36_sigSNP.hapmap -exportType Hapmap


#################################################################################
# Increasing BG SNPs

mkdir GAPIT_DArT_BGSNPs
cd GAPIT_DArT_BGSNPs

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --maf 0.1 --thin 5000 --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/iter_4_2_thin5000_GAPIT_DArT_BGSNPs # After filtering, kept 38823 out of a possible 33269031 Sites ; Run Time = 2202.00 seconds

#################################################################################
# GAPIT_39 - 3000 randsel all SNPs of putative genotype panel 
In R:
original_markers <- read.csv("~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/DArTag_form.csv") %>% select(Chrom, ChromPosPhysical, Comments)
selected_markers <- original_markers %>% group_by(Comments) %>% sample_n(size = 1000, replace = FALSE) %>% ungroup() %>% select(-c(Comments))
write.table(selected_markers, file = "~/RBGSyd_Technical Officer/MQuin/Outgroups/selected_markers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

In shell
cp "/home/karina/MQuin/Outgroups/data/meta/selected_markers.txt" /home/karina/mqgwas/iter_6/data/meta/

/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_4/gwas_4/Mq_filt_cat.beagle.vcf" --positions "/home/karina/mqgwas/iter_6/data/meta/selected_markers.txt" --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites

bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites.recode.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites.recode.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites.recode.vcf.gz
gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_39/MQ_iter4_genotyping3000-sites.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_39_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_39_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_39_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_39_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_39_sigSNP.hapmap -exportType Hapmap


#################################################################################
# GAPIT_40 - 1000 randsel MR SNPs of putative genotype panel 
In R:
markers_3000 <- read.table (file = "~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/selected_markers.txt")
dartag_form <- read.csv (file = "~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/DArTag_form.csv")
dartag_form_MR <- dartag_form %>% filter(Comments =="MR")
markers_3000_mr <- markers_3000[markers_3000$V2 %in% dartag_form_MR$ChromPosPhysical,]
write.table(markers_3000_mr, file = "~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/selected_markers_MR.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_4/gwas_4/Mq_filt_cat.beagle.vcf" --positions "/home/karina/mqgwas/iter_6/data/meta/selected_markers_MR.txt" --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites

bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites.recode.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites.recode.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites.recode.vcf.gz
gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_40/MQ_iter4_genotypingMR1000-sites.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_40_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_40_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_40_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_40_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_40_sigSNP.hapmap -exportType Hapmap

################################################################################
# GAPIT_41 - Prediction of 555 markers from GAPIT_39. Used to test against log10 precipitation from worldclim 0.5 degree.

Run GAPIT_climate.R

################################################################################
# GAPIT_42 - top 2000 MR related epistatic hits

mkdir /home/karina/mqgwas/iter_6/data/GAPIT_42
cp "/home/karina/mqgwas/iter_4/epistasis/whole_significant_GP/data/Mq_filt_cat_MRepistatic2000.recode.vcf" /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf

bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf.gz
rm /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf
gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_42/Mq_filt_cat_MRepistatic2000.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_42_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_42_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_41_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_42_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_42_sigSNP.hapmap -exportType Hapmap

################################################################################
# GAPIT_43 - top 2000 MR related epistatic hits
C:\Users\swirl\OneDrive\Documents\RBGSyd_Technical Officer\MQuin\Seedling GWAS\Filtering\Iteration 6\code\RGAPIT_43.R

################################################################################
# GAPIT_44 - BSLMM top 1000 SNPs; # GAPIT_45 - BSLMM top 500 SNPs 

mkdir GAPIT_44
cd GAPIT_44

head -n 1001 "/home/karina/mqgwas/iter_4/gwas_7/Output/Bayesian_Sparse_Linear_Mixed_Model/COI/COI_20240607_054609/highest_effects/COI_mod_sub_GWAS_phenotype_height_COItransformedfoursqrt.part1_Mq_filt_cat.beagle_top1%eff.csv" | awk -F',' 'NR > 1 {print $2 "\t" $5}' gwas_7_highesteff_1000SNPs.csv > gwas_7_highesteff_1000SNPs_pos.csv

head -n 501 "/home/karina/mqgwas/iter_4/gwas_7/Output/Bayesian_Sparse_Linear_Mixed_Model/COI/COI_20240607_054609/highest_effects/COI_mod_sub_GWAS_phenotype_height_COItransformedfoursqrt.part1_Mq_filt_cat.beagle_top1%eff.csv" | awk -F',' 'NR > 1 {print $2 "\t" $5}' > gwas_7_highesteff_500SNPs_pos.csv

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions gwas_7_highesteff_1000SNPs_pos.csv --recode --recode-INFO-all --out GAPIT_44_1000
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions gwas_7_highesteff_500SNPs_pos.csv --recode --recode-INFO-all --out GAPIT_45_500

bgzip -c GAPIT_44_1000.recode.vcf > GAPIT_44_1000.recode.vcf.gz
tabix -p vcf GAPIT_44_1000.recode.vcf.gz
rm GAPIT_44_1000.recode.vcf
gzip -d GAPIT_44_1000.recode.vcf.gz 

bgzip -c GAPIT_45_500.recode.vcf > GAPIT_45_500.recode.vcf.gz
tabix -p vcf GAPIT_45_500.recode.vcf.gz
rm GAPIT_45_500.recode.vcf
gzip -d GAPIT_45_500.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_44/GAPIT_44_1000.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_44_sigSNP.vcf
mv /home/karina/mqgwas/iter_6/data/GAPIT_45/GAPIT_45_500.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_45_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_44_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_44_sigSNP_sort.vcf" -fileType VCF
/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_45_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_45_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_44_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_44_sigSNP.hapmap -exportType Hapmap
/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_45_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_45_sigSNP.hapmap -exportType Hapmap

################################################################################
# GAPIT_48 - Final 3,300 SNPs for DArTag panel- seedlings; also ran with binary 0, 1, 2 option

mkdir GAPIT_48
cd GAPIT_48

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "positions_3300.txt" --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/GAPIT_48/MQ_final_iter4_DArTag3300-sites
## After filtering, kept 3005 out of a possible 33269031 Sites (vcf_pos.txt) - confirmed missing ~300 were climate SNPs + 1 MR SNP

## Numeric for GT
/usr/bin/vcftools --vcf MQ_final_iter4_DArTag3300-sites.recode.vcf --012 --out /home/karina/mqgwas/iter_6/data/GAPIT_48/MQ_final_iter4_DArTag3300_numeric

bgzip -c MQ_final_iter4_DArTag3300-sites.recode.vcf > MQ_final_iter4_DArTag3300-sites.recode.vcf.gz
tabix -p vcf MQ_final_iter4_DArTag3300-sites.recode.vcf.gz
rm MQ_final_iter4_DArTag3300-sites.recode.vcf
gzip -d MQ_final_iter4_DArTag3300-sites.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_48/MQ_final_iter4_DArTag3300-sites.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_48_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_48_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_48_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_48_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_48_sigSNP.hapmap -exportType Hapmap

## For MR and epi only, sites saved to DArTag_response_MR-epi.txt
MR_epi_snps <- read.table(file="~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/DArTag_response_MR-epi.txt")
myG_head <- myG[1,]
myG_filt <- myG %>% filter(V4 %in% MR_epi_snps$V2)
myG <- rbind(myG_head, myG_filt)

################################################################################
# GAPIT_49 - Final 3,300 SNPs for DArTag panel - parents

mkdir GAPIT_49
cd GAPIT_49

cp "/home/karina/mqgwas/iter_6/data/GAPIT_48/positions_3300.txt" .

/usr/bin/vcftools --gzvcf "/recer1/karina/parents/iter_1/gwas_5/Mq_filt_cat.beagle.vcf.gz" --positions "positions_3300.txt" --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/GAPIT_49/MQ_final_iter4_parents_DArTag3300-sites
## After filtering, kept 2890 out of a possible 47731087 Sites

bgzip -c MQ_final_iter4_parents_DArTag3300-sites.recode.vcf > MQ_final_iter4_parents_DArTag3300-sites.recode.vcf.gz
tabix -p vcf MQ_final_iter4_parents_DArTag3300-sites.recode.vcf.gz
rm MQ_final_iter4_parents_DArTag3300-sites.recode.vcf
gzip -d MQ_final_iter4_parents_DArTag3300-sites.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_49/MQ_final_iter4_parents_DArTag3300-sites.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_49_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_49_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_49_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_49_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_49_sigSNP.hapmap -exportType Hapmap

/usr/bin/vcftools --vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_49_sigSNP.vcf" --012 --out /home/karina/mqgwas/iter_6/data/GAPIT_49/MQ_final_iter4_parents_DArTag3300_num

/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/parent_iter_2/data/catvcf/Mq_filtmiss_cat.vcf.gz" --positions "positions_3300.txt" --012 --out /home/karina/mqgwas/iter_6/data/GAPIT_49/MQ_final_iter4_parents_DArTag3300-sites_nonimput
################################################################################
# GAPIT_50 - Rerunning with 1,000 SNPs with top associations P < 9.45E-06

mkdir GAPIT_50
cd GAPIT_50

# Cp SNPs 316 - 1316 from iter_4 gwas_2 to positions.txt
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --positions "positions.txt" --recode --recode-INFO-all --out MQ_iter4_midassoc

bgzip -c MQ_iter4_midassoc.recode.vcf > MQ_iter4_midassoc.recode.vcf.gz
tabix -p vcf MQ_iter4_midassoc.recode.vcf.gz
rm MQ_iter4_midassoc.recode.vcf
gzip -d MQ_iter4_midassoc.recode.vcf.gz 

mv /home/karina/mqgwas/iter_6/data/GAPIT_50/MQ_iter4_midassoc.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_50_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_50_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_50_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_50_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_50_sigSNP.hapmap -exportType Hapmap

################################################################################
# GAPIT_51 - What if we didn't know which SNPs were significant and just annotated the genome for disease proteins?

# Using all NLR & NBARC genes for prediction
# 1. Pull out all positions from the annotated genome - CDS
# 2. Pull out all SNPs from hapA
# 3. GAPIT

awk '$3 == "CDS" {print $1"\t"$4"\t"$5}' "/recer1/rgds/genomes/Mquin/ms_25082023/nbarc_hapA.gff3" > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nbarc_pos.txt
awk '$3 == "CDS" {print $1"\t"$4"\t"$5}' "/recer1/rgds/genomes/Mquin/ms_25082023/nlr_hapA.gff3" > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt

# Modifying chromosome names in NLR gff
sed 's/MQHAP1/MqA_/g' /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos_mod.txt
sed 's/\.01//g' /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos_mod.txt > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt
sed 's/\.02//g' /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos_mod.txt
sed 's/\.03//g' /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos_mod.txt > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt
rm /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos_mod.txt

cat /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nbarc_pos.txt /home/karina/mqgwas/iter_6/data/GAPIT_51/data/nlr_pos.txt > /home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos.bed

# Merge overlapping regions
echo -e "chrom\tchromStart\tchromStop" > "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos_merge.bed" # Add header
sort -k1,1 -k2,2n "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos.bed" > "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos_sort.bed"
/data/genomics/apps/bedtools/bedtools2/bin/mergeBed -i "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos_sort.bed" >> "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos_merge.bed"


# Pull out pos
/usr/bin/vcftools --gzvcf "/home/karina/mqgwas/iter_4/data/catvcf/Mq_filt_cat.vcf.gz" --bed "/home/karina/mqgwas/iter_6/data/GAPIT_51/data/all_pos_merge.bed" --recode --recode-INFO-all --out /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP #Read 3738 BED file entries. After filtering, kept 190473 out of a possible 33269031 Sites


bgzip -c /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP.recode.vcf > /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP.recode.vcf.gz
tabix -p vcf /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP.recode.vcf.gz
gzip -d /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP.recode.vcf.gz


mv /home/karina/mqgwas/iter_6/data/GAPIT_51/GAPIT_51_sigSNP.recode.vcf /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_51_sigSNP.vcf

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -SortGenotypeFilePlugin -inputFile /home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_51_sigSNP.vcf -outputFile "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_51_sigSNP_sort.vcf" -fileType VCF

/data/genomics/apps/tassel/run_pipeline.pl -Xms64G -Xmx64G -vcf "/home/karina/mqgwas/iter_6/data/sigSNP_vcf/GAPIT_51_sigSNP_sort.vcf" -export /home/karina/mqgwas/iter_6/data/hapmap/GAPIT_51_sigSNP.hapmap -exportType Hapmap

########################### NLR only

NLR_pos <- read.table("nlr_pos.txt")

remove(myGAPIT)
genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration,"_sigSNP.hapmap.hmp.txt")
myG <- read.csv(file = genetic_data, sep = "\t", header = FALSE)

# Initialize myG_NLR as a data frame with no rows but the same structure as myG
myG_NLR <- myG[1, ]  # Empty data frame with same columns as myG
chr_list <- unique(myG$V3[-1])

for (chr in chr_list) {
  myG_chrom <- myG %>% filter(V3 == chr)
  NLR_pos_chrom <- NLR_pos %>% filter(V1 == gsub("MQ", "Mq", chr))
  
  # Predefine a list to collect rows instead of using `rbind` in the loop
  selected_rows <- vector("list", nrow(myG_chrom))
  row_index <- 1
  
  for (i in seq_len(nrow(myG_chrom))) {
    if (i %% 100 == 0) {
      print(paste("Running", i, "out of", nrow(myG_chrom), "in", chr, "- SNPs included", (nrow(myG_NLR) - 1)))
    }
    
    pos <- as.numeric(myG_chrom$V4[i])
    
    # Vectorized filtering of NLR_pos_chrom based on pos +/- 10000 range
    NLR_pos_chrom_sub <- NLR_pos_chrom %>%
      filter(V3 < (pos + 10000), V2 > (pos - 10000))
    
    # Check if any row in NLR_pos_chrom_sub satisfies the positional condition
    sel <- any(pos > NLR_pos_chrom_sub$V2 & pos < NLR_pos_chrom_sub$V3)
    
    if (sel) {
      # Store selected row in list
      selected_rows[[row_index]] <- myG_chrom[i, ]
      row_index <- row_index + 1
    }
  }
  
  # Bind selected rows for this chromosome and append to myG_NLR
  selected_rows <- do.call(rbind, selected_rows[1:(row_index - 1)])  # Only bind rows we filled
  myG_NLR <- rbind(myG_NLR, selected_rows)
}

write.csv (myG_NLR, "myG_NLR.csv")

# Result in myG_NLR

myG_NLR$V3 <- chromosome_mapping[myG_NLR$V3] 
myG_NLR$V3[1] <- "chrom"
myGAPIT <- tryCatch(expr = {(GAPIT(Y=my_Y_train, G=myG_NLR, PCA.total=3, model=c("gBLUP")))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=my_Y_train, G=myG_NLR, PCA.total=3, model=c("gBLUP"), Geno.View.output = FALSE))})

# Extracting predictions + plotting
prediction=myGAPIT$Pred
unique(prediction$RefInf)

colnames(Y)[1] <- "Taxa"
prediction_gt <- left_join(prediction, Y)

prediction_gt_filtNA <- prediction_gt %>% filter(!is.na(COI))
PredvCOI_lm <- lm(prediction_gt_filtNA$COI ~ prediction_gt_filtNA$Prediction, na.action="na.exclude")
summary(PredvCOI_lm)
PredvCOI_resid <- resid(PredvCOI_lm)

p1 <- ggplot(prediction_gt, aes(x=COI, y=Prediction)) +
  geom_point(aes(colour = RefInf)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
p1

ggsave (p1, file = "PredictionvGT.jpg", height = 6.5, width = 6.5)

jpeg("Residuals.jpeg")
plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
dev.off()

# R2 values
lm <- lm(Prediction ~ COI, data = subset(prediction_gt))
R2 <- summary(lm)$adj.r.squared

lm <- lm(Prediction ~ COI, data = subset(prediction_gt, RefInf==2))
R2_testing <- summary(lm)$adj.r.squared

# Reliability
h2 = myGAPIT$h2
prediction_gt$Reliability = 1-prediction_gt$PEV/h2

mean(prediction_gt$Reliability)

p2 <- ggplot(prediction_gt, aes(y=Reliability, x=Prediction)) +
  geom_point(aes(colour = RefInf)) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (x = "Genomic Predicted COI^0.25", y = "Reliability") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
p2

iteration_reliability <- rbind(iteration_reliability, cbind(iteration, mean(prediction_gt$Reliability), R2, R2_testing))
ggsave (p2, file = "PredictionvPEV.jpg", height = 6.5, width = 6.5)
