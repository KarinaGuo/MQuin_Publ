
library(tidyverse)

############ Getting climate data
library(sp)
library(terra)
library(raster)

#clim <- geodata::worldclim_country("Australia", var="prec", res=0.5, path="~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/clim_data/")
clim <- raster("C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/clim_data/wc2.1_country/AUS_wc2.1_30s_prec.tif")

parent_meta_locs <- parent_meta %>% dplyr::select(longitude, latitude, NSWID) %>% filter(!is.na(latitude) & !is.na(NSWID))
coords <- parent_meta_locs %>% dplyr::select(longitude, latitude) %>% unique()
sp_points <- SpatialPoints(coords)
#sp_points_vect <- vect(sp_points, crs="+proj=longlat +datum=WGS84") # Comment out for downl load
prec <- extract(clim, sp_points)
#prec <- extract(clim, sp_points_vect)  # Comment out for downl load
coords_prec <- cbind(coords, prec) 
  #%>% dplyr::select(longitude,latitude,AUS_wc2.1_30s_prec_1)  # Comment out for downl load
colnames(coords_prec)[3] <- "precipitation"

parent_meta_climated <- left_join(parent_meta_locs, coords_prec)
parent_meta_climated <- parent_meta_climated %>% filter(!is.na(precipitation))  %>% unique() %>% dplyr::select(NSWID,precipitation)

# Appending to seedling ID
seedlings_all_subs <- seedlings_all_library %>% dplyr::select(LIBRARY, NSWID)
seedlings_all_subs_clim <- left_join(seedlings_all_library, parent_meta_climated) %>% unique()

## For parent locs
#parent_clim <- left_join(parent_meta_climated, parent_meta[,c("sample_lib_NSW", "NSWID")]) %>% dplyr::select(sample_lib_NSW, precipitation) %>% filter(!is.na(sample_lib_NSW)); colnames(parent_clim) <- c("LIBRARY", "precipitation")

################### Getting genotype data ready
# Climate SNPs from DArTag
DArTag_form <- read.csv(file="C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/DArTag/Melaleuca_quinquenervia_RBGSyd_DArTag_Marker_Form.csv")
DArTag_climate <- DArTag_form %>% filter(Comments=="climate") %>% dplyr::select(Chrom, ChromPosPhysical)


# Filtering hapmap
GT_data_raw <- read.csv(file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_39_sigSNP.hapmap.hmp.txt", sep = "\t", header = FALSE)
GT_data_raw_header <- GT_data_raw[1,]
GT_data_raw_climfilt <- GT_data_raw %>% filter(V4 %in% DArTag_climate$ChromPosPhysical)
GT_data_raw_climfilt <- rbind(GT_data_raw_header, GT_data_raw_climfilt)
myG <- GT_data_raw_climfilt

# climate_snps <- read.table(file="~/RBGSyd_Technical Officer/MQuin/Genotyping_Model/DArTag_response_climate.txt")
# GT_data_raw <- read.csv(file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_49_sigSNP.hapmap.hmp.txt", sep = "\t", header = FALSE)
# GT_data_raw_header <- GT_data_raw[1,]
# GT_data_raw_climfilt <- GT_data_raw %>% filter(V4 %in% climate_snps$V2)
# myG <- rbind(GT_data_raw_header, GT_data_raw_climfilt)

# Filtering climate seedlings to those with GT
seedlings_filt_clim <- seedlings_all_subs_clim %>% filter(LIBRARY %in% GT_data_raw[1,])
# parent_filt_clim <- parent_clim %>% filter(LIBRARY %in% GT_data_raw[1,])

# Prepping data
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

iteration=39
file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration)
if(!file.exists(file_dir)){
  print(paste("Directory", file_dir, "does not exist and is being created :)")) 
  dir.create(file_dir)
}
setwd(file_dir)

myG$V3 <- chromosome_mapping[myG$V3] 
myG$V3[1] <- "chrom"

Y <- seedlings_filt_clim %>% dplyr::select(LIBRARY, precipitation)
# Y <- parent_filt_clim %>% dplyr::select(LIBRARY, precipitation)
Y$precipitation <- log10(Y$precipitation)
my_Y_train <- Y[-c(1:50),]

source("http://zzlab.net/GAPIT/gapit_functions.txt")
myGAPIT <- tryCatch(expr = {(GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP")))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP"), Geno.View.output = FALSE))})

# Note only 555 SNPs

prediction=myGAPIT$Pred
unique(prediction$RefInf)

colnames(Y)[1] <- "Taxa"
prediction_gt <- left_join(prediction, Y)

prediction_gt_filtNA <- prediction_gt %>% filter(!is.na(precipitation))
PredvCOI_lm <- lm(prediction_gt_filtNA$precipitation ~ prediction_gt_filtNA$Prediction, na.action="na.exclude")
summary(PredvCOI_lm)
PredvCOI_resid <- resid(PredvCOI_lm)

jpeg("Residuals.jpeg")
plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
dev.off()

p1 <- ggplot(prediction_gt, aes(x=precipitation, y=Prediction)) +
  geom_point(aes(colour = RefInf)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (y = "Genomic Predicted log10(precipitation)", x = "Ground-truth log10(precipitation)") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
p1
ggsave (p1, file = "PredictionvGT.jpg", height = 6.5, width = 6.5)

lm <- lm(Prediction ~ precipitation, data = subset(prediction_gt))
R2 <- summary(lm)$adj.r.squared

lm <- lm(Prediction ~ precipitation, data = subset(prediction_gt, RefInf==2))
R2_testing <- summary(lm)$adj.r.squared

write.csv(prediction_gt, file="prediction_gt.csv", row.names = F)

#################################################################### Rerunning above for the final GT Dataset
