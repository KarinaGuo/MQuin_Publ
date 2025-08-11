convert_vcf012_gt <- function(dir, fn) {
  
  gt_fil <- paste0(dir,fn,".012")
  nm_fil <- paste0(dir,fn,".012.indv")
  po_fil <- paste0(dir,fn,".012.pos")
  
  # read
  gt <- read.table(gt_fil,header=F, sep="\t", row.names=1)
  gt[gt == -1] <- NA
  nm <- read.table(nm_fil,header=F, sep="\t")
  po <- read.table(po_fil,header=F, sep="\t")
  ppo <- paste(po[,1],po[,2],sep="_")
  
  rownames(gt) <- as.matrix(nm)
  colnames(gt) <- as.matrix(ppo)
  
  return(gt)
  
}


#################################### Set 1
setwd("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/ConStruct/Iter_1")

# Converting numeric thinned parent vcf to dart structure
dir = "data/"
fn = "Mq_vcf_thinnedSNPs_numeric"

gt <- as.matrix(convert_vcf012_gt(dir, fn))

# CLustering FID as populations
library(geosphere) 
library(cluster)

parent_meta_filtNA <- parent_meta %>% filter(!is.na(latitude))

distance_matrix <- as.matrix(cbind(parent_meta_filtNA$longitude,parent_meta_filtNA$latitude))
distance <- (distm(distance_matrix))
clusters <- as.hclust(agnes(distance, diss = T))
clust <- cutree(clusters, h = 1000)
distance_matrix <- as.data.frame(distance_matrix)
distance_matrix$cluster <- clust
colnames(distance_matrix) <- c("longitude", "latitude", "cluster")
distance_matrix_clustered <- distance_matrix %>% group_by (cluster) %>% summarise (mean_long = mean(longitude), mean_lat = mean(latitude))


# Manually creating a meta file for samples
parent_meta_filt <- left_join(parent_meta, distance_matrix, by = c('latitude', 'longitude'), multiple = "first") # clusters
parent_meta_filt <- parent_meta_filt %>% select(sample_lib_NSW, latitude, longitude, cluster)


parent_meta_filt <- parent_meta_filt[parent_meta_filt$sample_lib_NSW %in% rownames(gt),]
parent_meta_filt <- parent_meta_filt[order(parent_meta_filt$sample_lib_NSW),]
parent_meta_filt <- parent_meta_filt[!is.na(parent_meta_filt$latitude),]


analyses = NULL
RR1 = parent_meta_filt[!is.na(parent_meta_filt$latitude), "cluster"]
analyses$RR1 <- RR1


gt_filt <- gt[rownames(gt) %in% parent_meta_filt$sample_lib_NSW,] # Filter gt of those with no lat long info

parent_meta_filt_reorder <- parent_meta_filt[match(rownames(gt_filt), parent_meta_filt$sample_lib_NSW), ] # ensure matching order

ggplot() +
  geom_point (data = parent_meta_filt_reorder, aes(x=longitude, y=latitude, colour = as.character(cluster))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs (x = "Longitude", y = "Latitude") +
  theme_minimal ()

sample_names = parent_meta_filt_reorder$sample_lib_NSW

lat = parent_meta_filt_reorder$latitude
long = parent_meta_filt_reorder$long

site = paste(lat,long)

meta <- list(analyses, lat, long, sample_names, site)
names(meta) <- c("analyses", "lat", "long", "sample_names", "site")

data = list(gt_filt, sample_names, meta)
names(data) <- c("gt", "sample_names" ,"meta")

save(data, file = "data/VCFData_wmeta.Rdata")

# Generating population alleles
library(RRtools)
basedir = paste0(getwd(), "/")
species = "MelaleucaQuinquenervia"
dir.create("MelaleucaQuinquenervia")

treatment <- "run_1"
dataset = NA

# Below is from miras 'conStruct run.R' script
population_allele_stats  <- calculate.population.allele.stats(data, pop = data$meta$analyses$RR1)
population_spatial_dist  <- population.pw.spatial.dist(data, data$meta$analyses$RR1)

ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
if ( length(ind_NA_loci) > 0 ) {
  cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
  population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
  population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
  population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
  population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
}


dir <- paste(basedir, "output", sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen",sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen/",treatment,sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists...  \n")
}

cs_dir    <- paste(basedir,species,"/popgen/",treatment,"/conStruct", sep="")

if(!dir.exists(cs_dir)) {
  cat("  conStruc directory: ", cs_dir, " does not exist and is being created. \n")
  dir.create(cs_dir)
} else {
  cat("  conStruc directory: ", cs_dir, " already exists, content will be overwritten. \n")
}

cs_object_file   <- paste(cs_dir,"/",species,"_",dataset,".rda",sep="")

counts       <- population_allele_stats$minor
sample_sizes <- population_allele_stats$sample
lon_lat  <- population_spatial_dist$pop_info$lon_lat
freq <- counts/sample_sizes
prefix       <- paste(species, "_",dataset,sep="")

n <- nrow(lon_lat)
s <- mat.or.vec(n, n)

for (a in 1:n){
  lla <- c(lon_lat[a, 1], lon_lat[a, 2])
  for (b in 1:n){
    llb <- c(lon_lat[b, 1], lon_lat[b, 2])
    d <- distCosine(lla, llb)
    s[a, b] <- d
  }
}

cs <- list(counts=counts, sample_sizes=sample_sizes, lon_lat=lon_lat, freq=freq, dist=s, cs_dir=cs_dir, prefix=prefix)
save(cs, file=cs_object_file)

# Test
library(conStruct)
k=2

con <- conStruct(spatial= TRUE, K=2, freqs = cs$freq, geoDist = cs$dist, coords = cs$lon_lat, prefix = paste0("output/",cs$prefix,"_TRUE_K", k))


#################################### Set 2
# Set includes minim pop clustered as <6 ensure Mungo Brush Rd & Seal rocks Rd is one and Mile Lakes is independent. Add in QLD samples

# Converting numeric thinned parent vcf to dart structure
setwd("C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Seedling GWAS/ConStruct/Iter_2")
dir = "data/"
fn = "Mq_vcf_thinnedSNPs_numeric"

gt <- as.matrix(convert_vcf012_gt(dir, fn))

# Adding outgroup QLD
outgroups_QLD <- outgroups_merged %>% filter(species == "Melaleuca quinquenervia" & Site == "queensland") %>% filter(!is.na(sample_lib_NSW)) %>% select (sample_lib_NSW, latitude, longitude, pops)
colnames(outgroups_QLD) <- c("sample_lib_NSW", "latitude", "longitude", "Pop")

# Cleaning up the meta files
parent_meta_filt <- parent_meta %>% select("sample_lib_NSW", "latitude", "longitude", "Pop")
parent_meta_filt_out <- rbind(parent_meta_filt, outgroups_QLD)  %>% filter(!is.na(sample_lib_NSW))

parent_meta_filt_out %>% group_by(Pop) %>% summarise (n=n()) # One pop has 5, we'll live. Leaving qld ones as is

parent_meta_gtonly <- parent_meta_filt_out[parent_meta_filt_out$sample_lib_NSW %in% rownames(gt),]
parent_meta_reorder <- parent_meta_gtonly[order(parent_meta_gtonly$sample_lib_NSW),]
parent_meta_filtlat <- parent_meta_reorder[!is.na(parent_meta_reorder$latitude),]

# Manually creating a meta file for samples
gt_filt <- gt[rownames(gt) %in% parent_meta_filtlat$sample_lib_NSW,] # Filter gt of those with no lat long info

parent_meta_filt <- parent_meta_filtlat %>% filter(sample_lib_NSW %in% rownames(gt_filt))
parent_meta_filt_reorder <- parent_meta_filt[match(rownames(gt_filt), parent_meta_filt$sample_lib_NSW), ] # ensure matching order

analyses = NULL
RR1 = parent_meta_filt[!is.na(parent_meta_filt$latitude), "Pop"]
analyses$RR1 <- RR1




## Checking pops
ggplot() +
  geom_point (data = parent_meta_filt_reorder, aes(x=longitude, y=latitude, colour = as.character(Pop))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs (x = "Longitude", y = "Latitude") +
  theme_minimal ()


sample_names = parent_meta_filt_reorder$sample_lib_NSW

lat = parent_meta_filt_reorder$latitude
long = parent_meta_filt_reorder$long

site = paste(lat,long)

meta <- list(analyses, lat, long, sample_names, site)
names(meta) <- c("analyses", "lat", "long", "sample_names", "site")

data = list(gt_filt, sample_names, meta)
names(data) <- c("gt", "sample_names" ,"meta")

save(data, file = "data/VCFData_wmeta.Rdata")

# Generating population alleles
library(RRtools)
basedir = paste0(getwd(), "/")
species = "MelaleucaQuinquenervia"
dir.create("MelaleucaQuinquenervia")

treatment <- "run_2"
dataset = NA

# Below is from miras 'conStruct run.R' script
population_allele_stats  <- calculate.population.allele.stats(data, pop = data$meta$analyses$RR1)
population_spatial_dist  <- population.pw.spatial.dist(data, data$meta$analyses$RR1)

ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
if ( length(ind_NA_loci) > 0 ) {
  cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
  population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
  population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
  population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
  population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
}


dir <- paste(basedir, "output", sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen",sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen/",treatment,sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists...  \n")
}

cs_dir    <- paste(basedir,species,"/popgen/",treatment,"/conStruct", sep="")

if(!dir.exists(cs_dir)) {
  cat("  conStruc directory: ", cs_dir, " does not exist and is being created. \n")
  dir.create(cs_dir)
} else {
  cat("  conStruc directory: ", cs_dir, " already exists, content will be overwritten. \n")
}

cs_object_file   <- paste(cs_dir,"/",species,"_",dataset,".rda",sep="")

counts       <- population_allele_stats$minor
sample_sizes <- population_allele_stats$sample
lon_lat  <- population_spatial_dist$pop_info$lon_lat
freq <- counts/sample_sizes
prefix       <- paste(species, "_",dataset,sep="")

n <- nrow(lon_lat)
s <- mat.or.vec(n, n)

for (a in 1:n){
  lla <- c(lon_lat[a, 1], lon_lat[a, 2])
  for (b in 1:n){
    llb <- c(lon_lat[b, 1], lon_lat[b, 2])
    d <- distCosine(lla, llb)
    s[a, b] <- d
  }
}

cs <- list(counts=counts, sample_sizes=sample_sizes, lon_lat=lon_lat, freq=freq, dist=s, cs_dir=cs_dir, prefix=prefix)
save(cs, file=cs_object_file)


#################################### Set 3
# Update to 60k SNPs, remove QLD Pop 20645

# Converting numeric thinned parent vcf to dart structure
setwd("C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Seedling GWAS/ConStruct/Iter_3")
dir = "data/"
fn = "Mq_vcf_thinnedSNPs_numeric"

gt <- as.matrix(convert_vcf012_gt(dir, fn))

# Adding outgroup QLD
outgroups_QLD <- outgroups_merged %>% filter(species == "Melaleuca quinquenervia" & Site == "queensland") %>% filter(!is.na(sample_lib_NSW)) %>% select (sample_lib_NSW, latitude, longitude, pops) 
colnames(outgroups_QLD) <- c("sample_lib_NSW", "latitude", "longitude", "Pop")
outgroups_QLD <- outgroups_QLD %>% filter (Pop != "20645") # Remove QLD Pop 20645

# Cleaning up the meta files
parent_meta_filt <- parent_meta %>% select("sample_lib_NSW", "latitude", "longitude", "Pop")
parent_meta_filt_out <- rbind(parent_meta_filt, outgroups_QLD)  %>% filter(!is.na(sample_lib_NSW))

parent_meta_filt_out %>% group_by(Pop) %>% summarise (n=n()) # One pop has 5, we'll live. Leaving qld ones as is

parent_meta_gtonly <- parent_meta_filt_out[parent_meta_filt_out$sample_lib_NSW %in% rownames(gt),]
parent_meta_reorder <- parent_meta_gtonly[order(parent_meta_gtonly$sample_lib_NSW),]
parent_meta_filtlat <- parent_meta_reorder[!is.na(parent_meta_reorder$latitude),]

# Manually creating a meta file for samples
gt_filt <- gt[rownames(gt) %in% parent_meta_filtlat$sample_lib_NSW,] # Filter gt of those with no lat long info

parent_meta_filt <- parent_meta_filtlat %>% filter(sample_lib_NSW %in% rownames(gt_filt))
parent_meta_filt_reorder <- parent_meta_filt[match(rownames(gt_filt), parent_meta_filt$sample_lib_NSW), ] # ensure matching order

analyses = NULL
RR1 = parent_meta_filt[!is.na(parent_meta_filt$latitude), "Pop"]
analyses$RR1 <- RR1




## Checking pops
ggplot() +
  geom_point (data = parent_meta_filt_reorder, aes(x=longitude, y=latitude, colour = as.character(Pop))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs (x = "Longitude", y = "Latitude") +
  theme_minimal ()


sample_names = parent_meta_filt_reorder$sample_lib_NSW

lat = parent_meta_filt_reorder$latitude
long = parent_meta_filt_reorder$long

site = paste(lat,long)

meta <- list(analyses, lat, long, sample_names, site)
names(meta) <- c("analyses", "lat", "long", "sample_names", "site")

data = list(gt_filt, sample_names, meta)
names(data) <- c("gt", "sample_names" ,"meta")

save(data, file = "data/VCFData_wmeta.Rdata")

# Generating population alleles
library(RRtools)
basedir = paste0(getwd(), "/")
species = "MelaleucaQuinquenervia"
dir.create("MelaleucaQuinquenervia")

treatment <- "run_3"
dataset = NA

# Below is from miras 'conStruct run.R' script
population_allele_stats  <- calculate.population.allele.stats(data, pop = data$meta$analyses$RR1)
population_spatial_dist  <- population.pw.spatial.dist(data, data$meta$analyses$RR1)

ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
if ( length(ind_NA_loci) > 0 ) {
  cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
  population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
  population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
  population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
  population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
}


dir <- paste(basedir, "output", sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen",sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen/",treatment,sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists...  \n")
}

cs_dir    <- paste(basedir,species,"/popgen/",treatment,"/conStruct", sep="")

if(!dir.exists(cs_dir)) {
  cat("  conStruc directory: ", cs_dir, " does not exist and is being created. \n")
  dir.create(cs_dir)
} else {
  cat("  conStruc directory: ", cs_dir, " already exists, content will be overwritten. \n")
}

cs_object_file   <- paste(cs_dir,"/",species,"_",dataset,".rda",sep="")

counts       <- population_allele_stats$minor
sample_sizes <- population_allele_stats$sample
lon_lat  <- population_spatial_dist$pop_info$lon_lat
freq <- counts/sample_sizes
prefix       <- paste(species, "_",dataset,sep="")

n <- nrow(lon_lat)
s <- mat.or.vec(n, n)

for (a in 1:n){
  lla <- c(lon_lat[a, 1], lon_lat[a, 2])
  for (b in 1:n){
    llb <- c(lon_lat[b, 1], lon_lat[b, 2])
    d <- distCosine(lla, llb)
    s[a, b] <- d
  }
}

cs <- list(counts=counts, sample_sizes=sample_sizes, lon_lat=lon_lat, freq=freq, dist=s, cs_dir=cs_dir, prefix=prefix)
save(cs, file=cs_object_file)


#################################### Set 3
# Update to 60k SNPs, rerun using all NSW samples with no techreps and all 17  QLD indivs

setwd("C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Parent/ConStruct/Iter_4/")
dir = "data/"
fn = "Mq_60kthin_numeric"


gt <- as.matrix(convert_vcf012_gt(dir, fn))

# Adding outgroup QLD
outgroups_QLD <- outgroups_merged %>% filter(species == "Melaleuca quinquenervia" & Site == "queensland") %>% filter(!is.na(sample_lib_NSW)) %>% select (sample_lib_NSW, latitude, longitude, pops) 
colnames(outgroups_QLD) <- c("sample_lib_NSW", "latitude", "longitude", "Pop")

## Dropping samples with with less than 3 sampls
outgroups_QLD <- outgroups_QLD %>% filter (Pop != "20645") # Remove QLD Pop 20645

# Cleaning up the meta files
parent_meta_filt <- parent_meta %>% select("sample_lib_NSW", "latitude", "longitude", "Pop")
parent_meta_filt_out <- rbind(parent_meta_filt, outgroups_QLD)  %>% filter(!is.na(sample_lib_NSW))

parent_meta_filt_out %>% group_by(Pop) %>% summarise (n=n()) %>% arrange(n)# One pop has 5, we'll live. Leaving qld ones as is

parent_meta_gtonly <- parent_meta_filt_out[parent_meta_filt_out$sample_lib_NSW %in% rownames(gt),]
parent_meta_reorder <- parent_meta_gtonly[order(parent_meta_gtonly$sample_lib_NSW),]
parent_meta_filtlat <- parent_meta_reorder[!is.na(parent_meta_reorder$latitude),]

# Manually creating a meta file for samples
gt_filt <- gt[rownames(gt) %in% parent_meta_filtlat$sample_lib_NSW,] # Filter gt of those with no lat long info

parent_meta_filt <- parent_meta_filtlat %>% filter(sample_lib_NSW %in% rownames(gt_filt))
parent_meta_filt_reorder <- parent_meta_filt[match(rownames(gt_filt), parent_meta_filt$sample_lib_NSW), ] # ensure matching order

analyses = NULL
RR1 = parent_meta_filt[!is.na(parent_meta_filt$latitude), "Pop"]
analyses$RR1 <- RR1

## Checking pops
ggplot() +
  geom_point (data = parent_meta_filt_reorder, aes(x=longitude, y=latitude, colour = as.character(Pop))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs (x = "Longitude", y = "Latitude") +
  theme_minimal ()


sample_names = parent_meta_filt_reorder$sample_lib_NSW

lat = parent_meta_filt_reorder$latitude
long = parent_meta_filt_reorder$long

site = paste(lat,long)

meta <- list(analyses, lat, long, sample_names, site)
names(meta) <- c("analyses", "lat", "long", "sample_names", "site")

data = list(gt_filt, sample_names, meta)
names(data) <- c("gt", "sample_names" ,"meta")

save(data, file = "data/VCFData_wmeta.Rdata")

# Generating population alleles
library(RRtools)
basedir = paste0(getwd(), "/")
species = "MelaleucaQuinquenervia"
dir.create("MelaleucaQuinquenervia")


treatment <- "run_4"
dataset = NA

# Below is from miras 'conStruct run.R' script
population_allele_stats  <- calculate.population.allele.stats(data, pop = data$meta$analyses$RR1)
population_spatial_dist  <- population.pw.spatial.dist(data, data$meta$analyses$RR1)

ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
if ( length(ind_NA_loci) > 0 ) {
  cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
  population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
  population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
  population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
  population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
}


dir <- paste(basedir, "output", sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen",sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen/",treatment,sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists...  \n")
}

cs_dir    <- paste(basedir,species,"/popgen/",treatment,"/conStruct", sep="")

if(!dir.exists(cs_dir)) {
  cat("  conStruc directory: ", cs_dir, " does not exist and is being created. \n")
  dir.create(cs_dir)
} else {
  cat("  conStruc directory: ", cs_dir, " already exists, content will be overwritten. \n")
}

cs_object_file   <- paste(cs_dir,"/",species,"_",dataset,".rda",sep="")

counts       <- population_allele_stats$minor
sample_sizes <- population_allele_stats$sample
lon_lat  <- population_spatial_dist$pop_info$lon_lat
freq <- counts/sample_sizes
prefix       <- paste(species, "_",dataset,sep="")

n <- nrow(lon_lat)
s <- mat.or.vec(n, n)

for (a in 1:n){
  lla <- c(lon_lat[a, 1], lon_lat[a, 2])
  for (b in 1:n){
    llb <- c(lon_lat[b, 1], lon_lat[b, 2])
    d <- distCosine(lla, llb)
    s[a, b] <- d
  }
}

cs <- list(counts=counts, sample_sizes=sample_sizes, lon_lat=lon_lat, freq=freq, dist=s, cs_dir=cs_dir, prefix=prefix)
save(cs, file=cs_object_file)

