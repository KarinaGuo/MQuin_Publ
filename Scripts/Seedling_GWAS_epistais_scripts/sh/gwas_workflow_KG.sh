  BASE_DIR="/home/karina/mqgwas/iter_4"

# Setting up directory structure and sym-linking relevant raw files do local directory 

  mkdir -p ${BASE_DIR}/data/meta/ 
  mkdir -p ${BASE_DIR}/data/rawvcf/ 
  mkdir -p ${BASE_DIR}/data/filtvcf_DP6_GQ20/
  mkdir -p ${BASE_DIR}/data/filtvcf_indels_biallelic/
  mkdir -p ${BASE_DIR}/data/filtvcf_exlowdepth/
  mkdir -p ${BASE_DIR}/data/filtvcf_missing/
  mkdir -p ${BASE_DIR}/data/qualvcf/
  mkdir -p ${BASE_DIR}/data/catvcf/
  mkdir -p ${BASE_DIR}/data/filtvcf_hwe/
  
  cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 96 ${BASE_DIR}/code/sh/submit_mpileup_script.sh {} "${BASE_DIR}"

# Removing indels and only keeping biallelic sites
  cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 96 ${BASE_DIR}/code/sh/submit_vcftools_script.sh {} "${BASE_DIR}"

# The list of individuals to be filtered per depth has already been calculated previously therefore left out from recalculating


# Remove all low depth sites. Masking sites with DP6 GQ20. Removing by missingness 200.
  cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 96 ${BASE_DIR}/code/sh/submit_vcftools_script_2.sh {} "${BASE_DIR}"

  cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 96 ${BASE_DIR}/code/sh/submit_vcftools_hwe.sh {} "${BASE_DIR}"

/data/genomics/apps/bcftools-1.17/bcftools concat -f ${BASE_DIR}/data/meta/Mqui_region_vcfs_hwe.txt -Oz -o ${BASE_DIR}/data/catvcf/Mq_filthwe_cat.vcf.gz

java -jar -Xmx900g /data/genomics/apps/beagle/beagle.22Jul22.46e.jar gt=${BASE_DIR}/data/catvcf/Mq_filthwe_cat.vcf.gz out=${BASE_DIR}/data/catvcf/Mq_hwe_cat.beagle impute=true nthreads=96

# Unzip imputed file and moved relevant GWAS file
gzip -d "${BASE_DIR}/data/catvcf/Mq_hwe_cat.beagle.vcf.gz"
mv ${BASE_DIR}/data/catvcf/Mq_hwe_cat.beagle.vcf ${BASE_DIR}/gwas

# Running GWAS on COI and RGR!
docker run -v ${BASE_DIR}/gwas/:/vcf2gwas/ fvogt257/vcf2gwas -v Mq_hwe_cat.beagle.vcf -pf GWAS_phenotype_height_COItransformedfoursqrt.csv -ap -lmm --minaf 0.05








######## OLD
# Adding relevant files. Remember to move in the following files and change file paths of 2): 
#1) libarary_list.txt to /data/meta (the list of individuls phenotyped). 
#2) Mqui_v1_hap_regions.tsv to meta (list of rawvcf regions in order)
#3) Mqui_region_vcfs.txt to meta (list of regions in order after filtering for missing, in the following directory${BASE_DIR}/data/filtvcf_missing/). 
#4) GWAS_phenotype_height_COI.csv to gwas (phenotype list by library ID)

ls -d /home/jgb/mqgwas/data/rawvcf/* > ${BASE_DIR}/rawvcf_list.txt
for files in $(cat ${BASE_DIR}/rawvcf_list.txt); do ln -s $files ${BASE_DIR}/data/rawvcf/; done

### Filtering

# Masking sites with DP6 GQ20, removing individuals not phenotyped, removing indels and only keeping biallelic sites, and calculating mean depth of individuals
cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 48 ${BASE_DIR}/code/sh/submit_vcftools_script.sh {} "${BASE_DIR}"

# From the output of mean depth per individual, concatenate the metric output file into one file
mv "${BASE_DIR}/data/qualvcf/mq.depth.MqA_CHR01:1-5000000.out.idepth" "${BASE_DIR}/data/evalmetric_Mq_depth_regiondepth_file.out.idepth"

ls -d ${BASE_DIR}/data/qualvcf/* > ${BASE_DIR}/data/meta/depth_file_list.txt

for file in $(cat ${BASE_DIR}/data/meta/depth_file_list.txt); do
    tail -n +2 "$file" >> "${BASE_DIR}/data/evalmetric_Mq_depth_phenosub_file.out.idepth"
done

# In R filter the depth metrics using process_depth_data.R which is:
Rscript --vanilla ${BASE_DIR}/code/r/process_depth_data.R  ${BASE_DIR}

# Removed all low depth sites. Removing by missingness 150 & recalculated depth
cat ${BASE_DIR}/data/meta/Mqui_v1_hap_regions.tsv | parallel -j 48 ${BASE_DIR}/code/sh/submit_vcftools_script_2.sh {} "${BASE_DIR}"

### Imputing

#Concatenated vcf files filtered for depth and missingness
#Copied over /home/jgb/mqgwas/data/meta/Mqui_region_vcfs.txt which was generated using /home/jgb/mqgwas/code/py/region2vcffilenames.py and changed file paths

/data/genomics/apps/bcftools-1.17/bcftools concat -f ${BASE_DIR}/data/meta/Mqui_region_vcfs.txt -Oz -o ${BASE_DIR}/data/catvcf/Mq_filtmiss_cat.vcf.gz

# beagle: impute sporadic missing genotypes
java -jar -Xmx900g /data/genomics/apps/beagle/beagle.22Jul22.46e.jar gt=${BASE_DIR}/data/catvcf/Mq_filtmiss_cat.vcf.gz out=${BASE_DIR}/data/catvcf/Mq_filt_cat.beagle impute=true nthreads=48

###
# Unzip imputed file and moved relevant GWAS file
gzip -d "${BASE_DIR}/data/catvcf/Mq_filt_cat.beagle.vcf.gz"
mv ${BASE_DIR}/data/catvcf/Mq_filt_cat.beagle.vcf.gz ${BASE_DIR}/gwas

# Running GWAS on COI and RGR!
docker run -v ${BASE_DIR}/gwas/:/vcf2gwas/ fvogt257/vcf2gwas -v Mq_filt_cat.beagle.vcf -pf GWAS_phenotype_height_COI.csv -ap -lmm
