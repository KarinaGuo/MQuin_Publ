library(tidyverse)
library(geosphere) 
library(cluster)

# Clustering FID by population
# Removing FID not in seedling dataset
parent_meta_sub <- subset(parent_meta, parent_meta$FID %in% COI_height_IID$FID)
parent_meta_sub <- parent_meta_sub %>% filter(!is.na(latitude))

# Calculating clusters
parent_meta_sub_matrix <- as.matrix(cbind(parent_meta_sub$longitude,parent_meta_sub$latitude))
distance <- (distm(parent_meta_sub_matrix))
clusters <- as.hclust(agnes(distance, diss = T))
clust <- cutree(clusters, h = 3000)

# Distinguishing each cluster 
parent_meta_sub$clust <- clust

# Plotting to check
tmp_seq <- seq(1, max(clust))
random_colors <- rgb(runif(max(clust)), runif(max(clust)), runif(max(clust)))
tmp_colours <- as.data.frame(cbind(tmp_seq, random_colors))
colnames(tmp_colours)[1] <- "clust"
tmp_colours$clust <- as.numeric(tmp_colours$clust)
parent_meta_sub_cols <- left_join(parent_meta_sub, tmp_colours)

ggplot(data= parent_meta_sub_cols, aes(x=longitude, y=latitude, col = random_colors)) + 
  geom_point() +
  theme_minimal() 

# Joining to seedling data
COI_height_IID_clustFID_coloured <- left_join (COI_height_IID, parent_meta_sub_cols, by = "FID") %>% select (LIBRARY, FID, latitude, longitude, clust, random_colors) %>% filter (!is.na(LIBRARY))

ggplot(data= COI_height_IID_clustFID_coloured, aes(x=longitude, y=latitude, color = random_colors, label = clust)) + 
  geom_text() +
  theme_minimal()

samples_testing <- subset(COI_height_IID_clustFID_coloured, COI_height_IID_clustFID_coloured$clust %in% c(3, 2, 6, 9, 8))
# 152 samples

# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
library(GAPIT)
library(tidyverse)

#Import files
myY <- mq_phenotype_height_COI
myY$COI <- (myY$COI)^0.25
Y=myY[,1:2]  
my_Y_train <- subset(Y, !(Y$LIBRARY %in% samples_testing$LIBRARY))
# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)


iteration_reliability <- data.frame()


  remove(myGAPIT)
  file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",25)
  if(!file.exists(file_dir)){
    print(paste("Directory", file_dir, "does not exist and is being created :)")) 
    dir.create(file_dir)
  }
  wd = file_dir
  genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",17,"_sigSNP.hapmap.hmp.txt")
  
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
  
  ggsave (p1, file = "PredictionvGT.jpg", height = 6.5, width = 6.5)
  
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
  
  iteration_reliability <- rbind(iteration_reliability, cbind(25, mean(prediction_gt$Reliability), R2))
  ggsave (p2, file = "PredictionvPEV.jpg", height = 6.5, width = 6.5)


prediction_gt$Prediction[prediction_gt$Prediction < 0] <- 0


iteration_reliability







