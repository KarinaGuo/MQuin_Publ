### Training North, Predicting South and viceversa


# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#library(GAPIT)
library(tidyverse)

#Import files
source("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Calculating metrics.R")

COI_height_IID_loc <- left_join(COI_height_IID %>% filter(!is.na(LIBRARY)), full_table_trim[,c("FID", "latitude", "longitude")], by = "FID")
COI_height_IID_loc$COI <- (COI_height_IID_loc$COI)^0.25

# Selected at natural separation points to balance minimal training data -29oS was selected for Northern training, -33oS was selected for Southern training

myY_south <- COI_height_IID_loc %>% 
  filter (!is.na(latitude) & latitude > -33) %>% 
  select(LIBRARY, COI) # Taking all those with a latitude higher than 33 to test the southern most population

myY_north <- COI_height_IID_loc %>% 
  filter (!is.na(latitude) & latitude < -29) %>% 
  select(LIBRARY, COI) # Taking all those with a latitude lower than 29 to test the northern most population

plot(COI_height_IID_loc[,c("longitude", "latitude")])
abline(h=-29)
abline(h=-33)
dev.off()

Y <- COI_height_IID_loc %>% select(LIBRARY, COI)

#my_Y_train <- myY_north
my_Y_train <- myY_south

# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

iterations <- c(48)

# For final GP selection in iter 48
MR_epi_snps <- read.table(file="~/RBGSyd_Technical Officer/MQuin/DArTag/DArTag_response_MR-epi.txt")

genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration,"_sigSNP.hapmap.hmp.txt")
myG <- read.csv(file = genetic_data, sep = "\t", header = FALSE)
myG$V3 <- chromosome_mapping[myG$V3] 
myG$V3[1] <- "chrom"

myG_head <- myG[1,]
myG_filt <- myG %>% filter(V4 %in% MR_epi_snps$V2)
myG <- rbind(myG_head, myG_filt)

iteration_reliability <- data.frame()

for (iteration in iterations){
  # Setting file paths
  remove(myGAPIT)
  file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration)
  if(!file.exists(file_dir)){
    print(paste("Directory", file_dir, "does not exist and is being created :)")) 
    dir.create(file_dir)
  }
  wd = file_dir
  
  # Running gBLUP for GP
  setwd(wd)
  myGAPIT <- tryCatch(expr = {GAPIT(Y = my_Y_train, G = myG, PCA.total = 3, model = c("gBLUP"), file.output = FALSE)}, error = function(e) {print("Running Geno.View.output = FALSE"); GAPIT(Y = my_Y_train, G = myG, PCA.total = 3, model = c("gBLUP"), Geno.View.output = FALSE)})
  
  # Extracting predictions + plotting
  prediction=myGAPIT$Pred
  table(prediction$RefInf)
  
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
  
  # ggsave (p1, file = "PredictionvGT.jpg", height = 6.5, width = 6.5)
  # 
  # jpeg("Residuals.jpeg")
  # plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
  # dev.off()
  
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
  #ggsave (p2, file = "PredictionvPEV.jpg", height = 6.5, width = 6.5)
}


prediction_gt$Prediction[prediction_gt$Prediction < 0] <- 0


iteration_reliability


p3 <- ggplot(subset(prediction_gt, RefInf==2), aes(x=COI, y=Prediction)) +
  geom_point(aes(colour = RefInf)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
  theme_bw() +
  theme(legend.position = "right") 
#scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
#ggsave (p3, file = "PredictionvGT_testing.jpg", height = 6.5, width = 6.5)

p1 + p2 + p3 + plot_layout(guides = "collect")
