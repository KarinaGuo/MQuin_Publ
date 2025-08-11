# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
library(GAPIT)
library(tidyverse)

#Import files
myY <- mq_phenotype_height_COI
myY$COI <- (myY$COI)^0.25
Y=myY[,1:2]  
tmpY_binaryCOI <- ifelse(Y$COI<1, 0, 1)
tmpY <- data.frame(cbind(myY[,1], tmpY_binaryCOI))
colnames(tmpY) <- c('Taxa', 'binary_COI')
tmpY <- subset(tmpY, !is.na(binary_COI))

tmpY_0 <- (subset(tmpY, binary_COI==0))[-c(1:25),]
tmpY_1 <- (subset(tmpY, binary_COI==1))[-c(1:25),]

my_Y_train <- rbind(tmpY_0, tmpY_1)
# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

iterations <- c(33)

  # Setting file paths
  remove(myGAPIT)
  file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_test")
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
  myGAPIT <- tryCatch(expr = {(GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP")))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP"), file.output=F))})
  
  # Extracting predictions + plotting
  prediction=myGAPIT$Pred
  unique(prediction$RefInf)
  
  colnames(Y)[1] <- "Taxa"
  prediction_gt <- left_join(prediction, tmpY)
  
  prediction_gt_filtNA <- prediction_gt %>% filter(!is.na(COI))
  PredvCOI_lm <- lm(prediction_gt_filtNA$COI ~ prediction_gt_filtNA$Prediction, na.action="na.exclude")
  summary(PredvCOI_lm)
  PredvCOI_resid <- resid(PredvCOI_lm)
  
  p1 <- ggplot() +
    geom_boxplot(data = prediction_gt, aes(x=binary_COI, y=Prediction)) +
    geom_jitter(data = prediction_gt, aes(x=binary_COI, y=Prediction, colour = RefInf), width = 0.1, alpha = 0.5) +
    labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))+
    labs(title="Training & testing")
  p1 
  
  p2 <- ggplot() +
    #geom_boxplot(data = subset(prediction_gt, RefInf==2), aes(x=binary_COI, y=Prediction)) +
    geom_point(data = prediction_gt, aes(x=0, y=Prediction, colour = binary_COI), alpha = 0.5) +
    labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
    theme_bw() +
    theme(legend.position = "bottom") 
    #scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing")) 
  p2
  
  p3 <- ggplot() +
    geom_boxplot(data = subset(prediction_gt, RefInf==2), aes(x=binary_COI, y=Prediction)) +
    geom_jitter(data = subset(prediction_gt, RefInf==2), aes(x=binary_COI, y=Prediction, colour = RefInf), width = 0.1, alpha = 0.5) +
    labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing")) +
    labs(title="Only testing")
  p3
  
  p4 <- ggplot() +
    #geom_boxplot(data = subset(prediction_gt, RefInf==2), aes(x=binary_COI, y=Prediction)) +
    geom_point(data = subset(prediction_gt, RefInf==2), aes(x=0, y=Prediction, colour = binary_COI), alpha = 0.5) +
    labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
    theme_bw() +
    theme(legend.position = "none") 
  #scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
  p4  

  library(patchwork)
  
  (p1+p2)/(p3+p4)
  
