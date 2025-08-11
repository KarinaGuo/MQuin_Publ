### Annotating results with resistance of Pegg

# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
library(GAPIT)
library(tidyverse)

#Import files
myY <- mq_phenotype_height_COI
myY$COI <- (myY$COI)^0.25
Y=myY[,1:2]  
my_Y_train <- Y[]
# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

iteration <- c(23)

  # Setting file paths
remove(myGAPIT)
file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_Pegg")
if(!file.exists(file_dir)){
  print(paste("Directory", file_dir, "does not exist and is being created :)")) 
  dir.create(file_dir)
}
wd = file_dir
genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration,"_sigSNP.hapmap.hmp.txt")

  # Running gBLUP for GP
setwd(wd)
myG <- read.csv(file = genetic_data, sep = "\t", header = FALSE)
myG$V3 <- chromosome_mapping[myG$V3] 
myG$V3[1] <- "chrom"
myGAPIT <- tryCatch(expr = {(GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP")))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP"), Geno.View.output = FALSE))})
  
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
  
jpeg("Residuals.jpeg")
plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
dev.off()
  
  # R2 values
lm <- lm(Prediction ~ COI, data = subset(prediction_gt))
R2 <- summary(lm)$adj.r.squared
  
  # Reliability
h2 = myGAPIT$h2
prediction_gt$Reliability = prediction_gt$PEV/h2
  
mean(prediction_gt$Reliability)
  
p2 <- ggplot(prediction_gt, aes(y=Reliability, x=Prediction)) +
    geom_point(aes(colour = RefInf)) +
    geom_smooth (method = "lm", se = F, linewidth = 0.5) +
    labs (x = "Genomic Predicted COI^0.25", y = "Reliability") +
    theme_bw() +
    theme(legend.position = "right") +
    scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
p2
  
iteration_reliability <- rbind(iteration_reliability, cbind(iteration, mean(prediction_gt$Reliability), R2))


############
#Annotating prediction plot

# Datasets of interest
# parent_mquin_rust  
# prediction + (sample_lib_full_outtrimmed + outgroups_merged)

prediction_inf <- prediction %>% filter(RefInf==2)

colnames(sample_lib_full_outtrimmed)[3] <- "Taxa"
prediction_FID <- left_join(prediction_inf, sample_lib_full_outtrimmed)
prediction_loc <- left_join(prediction_FID, outgroups_merged)
prediction_loc$latitude <- round(prediction_loc$latitude, 2)
prediction_loc$longitude <- round(prediction_loc$longitude,2)

Pegg_data <- read.table (file = "~/RBGSyd_Technical Officer/MQuin/Processing Meta/Pegg2018.txt", header = T)
colnames(Pegg_data) <- c("latitude", "longitude", "Pegg_SeedlingResistance")
Pegg_data$latitude <- round(Pegg_data$latitude,2)
Pegg_data$longitude <- round(Pegg_data$longitude,2)


prediction_loc_Pegg <- left_join(prediction_loc, Pegg_data)

ggplot() +
  geom_boxplot(data = prediction_loc_Pegg, aes(y=Prediction)) +
  geom_point(data = prediction_loc_Pegg, aes(x=0, y=Prediction)) +
  ggrepel::geom_label_repel(data = prediction_loc_Pegg, aes(x=0, y=Prediction, label = Pegg_SeedlingResistance, color = Pegg_SeedlingResistance), nudge_x = 0.15, nudge_y = 0.1, min.segment.length = 0.01) +
  scale_color_continuous(type = "viridis", name = "Pegg Resistant Seedling (%)") +
  theme_bw ()







