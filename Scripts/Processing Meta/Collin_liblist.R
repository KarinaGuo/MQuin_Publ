library(tidyverse)

# Reading in datasets
full_table <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Full_table.csv") # Maternal ID and the number of seedlings (n) grown 
colnames(full_table) <- c("FID", "SeedlingNumb", "species", "latitude", "longitude", "Pop", "ID", "IID")
outgroups <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Outgroups.csv") #Manually curated outgroup list from other datasets and external info from Jason G Bragg
colnames(outgroups) <- c("SampleNumb", "WellID", "SampleName", "NSWID", "Concentration", "species", "Site", "pop")
sample_list_mquin_PBI <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/sample_list_mquin_PBI.csv") # IID is the defining ID for each seedling. Heights (cm) were measured at two points at 23/08/23 and 16/09/23. FID is the maternal ID.
parent_mquin_rust <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/mquin_rust_phenotypic_data_20201218.csv") # Rust data of parent individuals
colnames(parent_mquin_rust) <- c("species", "NSWID", "FID", "latitude", "longitude", "height", "rust_infection", "LeavesShootInfect", "LeavesShootRustSchore1.5", "MinorBranchDB", "MajorBranchDB", "CrownDensity", "Notes")
parent_mquin_rust$LeavesShootInfect <- gsub("1-May", "1-5", parent_mquin_rust$LeavesShootInfect)
parent_mquin_rust$LeavesShootInfect <- gsub("Oct-50", "10-50", parent_mquin_rust$LeavesShootInfect)
parent_mquin_rust$LeavesShootInfect <- gsub("5-Oct", "5-10", parent_mquin_rust$LeavesShootInfect)
seedling_mquin_rust <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/mq_phenotypes.csv") # Rust data of seedling individuals
sample_lib <- read.table("~/RBGSyd_Technical Officer/MQuin/Processing Meta/geno_library_list.txt", header=TRUE) # sample_lib has the library names of all samples that were genotyped. They are named in the format of S_IID_Sample# or NSWID_IID_Sample#
plate_1 <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/plates 1-7 library prep(1).csv") # 1. Meta data for the library names, includes sample number and sample ID used to generate the library names
plate_2 <- read.csv("~/RBGSyd_Technical Officer/MQuin/Processing Meta/plate 8 library prep(1).csv")# 2. Meta data for the library names, includes sample number and sample ID used to generate the library names
seedlings_all <- read.csv ("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Mquin_samples_pheno.csv")# Meta data for seedlings including height and COI of all seedlings, including those that didn't get genotypes
colnames(seedlings_all) <- c("IID", "HT1", "HT2", "COI", "NSWID", "n","species", "latitude", 'longitude', "pop", "ID")


### Creating master datasets
## Parent datasets: NSWID, Library name, species name, infection, library details, outgroups

# Trimming full_table to one per maternal individual
full_table_trim <- full_table %>% group_by(FID) %>% filter(IID == max(IID)) %>% ungroup()
# Merging sample_lib and full_table_trim by FID first working with libraries of NSW_* then NSW*
sample_lib_NSW_ <- grep("NSW_", sample_lib$LIBRARY, value = TRUE)
FID <- sub("^NSW_(\\d+)_.*", "\\1", sample_lib_NSW_)
NSWID <- paste0("NSW", FID)
FID <- paste0("NSW_", FID)
sample_lib_NSW_ID <- as.data.frame(cbind(NSWID, FID, sample_lib_NSW_))
colnames(sample_lib_NSW_ID) <- c("NSWID", "FID", "sample_lib_NSW")

sample_lib_NSW <- grep("^NSW\\d+[^_]", sample_lib$LIBRARY, value = TRUE)
FID <- sub("^NSW(\\d+).*", "\\1", sample_lib_NSW)
NSWID <- paste0("NSW", FID)
FID <- paste0("NSW_", FID)
sample_lib_NSWID <- as.data.frame(cbind(NSWID, FID, sample_lib_NSW))
sample_lib_NSWBoth <- rbind(sample_lib_NSW_ID, sample_lib_NSWID)
sample_lib_full <- merge(sample_lib_NSWBoth, full_table_trim, by="FID", all=TRUE)
parent_meta <- merge(sample_lib_full, parent_mquin_rust, by=c("FID", "species", "longitude", "latitude"), all=TRUE)

# Outgroup dataset
outgroups_NSW <- grep("^NSW\\d+[^_]", outgroups$NSWID, value = TRUE)
FID <- sub("^NSW(\\d+).*", "\\1", outgroups_NSW)
NSWID <- paste0("NSW", FID)
FID <- paste0("NSW_", FID)
outgroups_NSWID <- as.data.frame(cbind(NSWID, FID))
sample_lib_full_outtrimmed <- sample_lib_full %>% select (FID, NSWID, sample_lib_NSW)
outgroups_geno <- merge(outgroups_NSWID, sample_lib_full_outtrimmed, by = c("FID", "NSWID"))
outgroups_merged <- left_join(outgroups_geno, outgroups, by = "NSWID") %>% mutate(source="outgroups")
outgroups_merged_trim <- outgroups_merged %>% select(sample_lib_NSW, species) %>% 
  filter (species == "Melaleuca quinquenervia")

outgroups_merged_trim[["species"]] <- "Melaleuca quinquenervia_QLD"

# NSW Parent datasets
trim_parent_meta <-  parent_meta %>% 
  filter(!is.na(sample_lib_NSW)) %>% 
  select(sample_lib_NSW, species)
trim_parent_meta <- trim_parent_meta[!(trim_parent_meta$sample_lib_NSW %in% outgroups_merged_trim$sample_lib_NSW), ]

# Lib name output list
merged_out_list <- rbind(trim_parent_meta, outgroups_merged_trim)

write.table(merged_out_list$sample_lib_NSW, file= "~/RBGSyd_Technical Officer/MQuin/Processing Meta/Collin_geno_out_library_list.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
