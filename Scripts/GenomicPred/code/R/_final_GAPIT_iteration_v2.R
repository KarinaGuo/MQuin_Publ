# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#library(GAPIT)
library(tidyverse)

# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

myY <- seedlings_all_library
myY$COI <- (myY$COI)^0.25
Y=myY[,c("LIBRARY","COI", "pop", "NSWID")]  
iterations <- c(48) # only 1 please

GAPIT_function <- function(iteration){
  remove(myGAPIT)
  remove(iteration_reliability)
  iteration_reliability <- data.frame()
  
  file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration)
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
  lm <- lm(Prediction ~ COI, data = subset(prediction_gt))
  R2 <- summary(lm)$adj.r.squared
  
  lm <- lm(Prediction ~ COI, data = subset(prediction_gt, RefInf==2))
  R2_testing <- summary(lm)$adj.r.squared
  
  # Reliability
  h2 = myGAPIT$h2
  prediction_gt$Reliability = 1-prediction_gt$PEV/h2
  
  mean(prediction_gt$Reliability)
  
  iteration_reliability <- rbind(iteration_reliability, cbind(iteration, mean(prediction_gt$Reliability), R2, R2_testing))
  
  GAPIT_run_res <- list(iteration_reliability, myGAPIT)
  
  return(GAPIT_run_res)
} 
  
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
    
  p1.5 <- ggplot(prediction_gt[prediction_gt$RefInf==2,], aes(x=COI, y=Prediction)) +
    geom_point(aes(colour = RefInf)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth (method = "lm", se = F, linewidth = 0.5) +
    labs (y = "Genomic Predicted COI^0.25", x = "Ground-truth COI^0.25") +
    theme_bw() +
    theme(legend.position = "right") +
    scale_color_discrete(name = "Prediction Category", labels = c("Testing"))
  p1.5
  
  ggsave (p1.5, file = "PredictionvGT_training.jpg", height = 6.5, width = 6.5)
  
  jpeg("Residuals.jpeg")
  plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
  dev.off()
  
  p2 <- ggplot(prediction_gt, aes(y=Reliability, x=Prediction)) +
    geom_point(aes(colour = RefInf)) +
    geom_smooth (method = "lm", se = F, linewidth = 0.5) +
    labs (x = "Genomic Predicted COI^0.25", y = "Reliability") +
    theme_bw() +
    theme(legend.position = "right") +
    scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
  p2

#Phenotype files
# 1. Validating by random pops
pop_all <- unique(seedlings_all_library$pop)
runs = 5
iteration_res_raw = NULL

for (iteration in iterations) {
  for (i in 1:runs){
    selected_populations <- sample(pop_all, 3)

    print(paste("Running pop:", list(selected_populations), ". On run number", i, "out of", runs)) 
    
    my_Y_train <- Y[!(Y$pop %in% selected_populations),1:2]
    iteration_pop_res = NULL
    
    GAPIT_run_res <- GAPIT_function(iteration)
    
    myGAPIT <- GAPIT_run_res[[2]]
    iteration_reliability <- GAPIT_run_res[[1]]
    
    iteration_pop_res <- as.data.frame(cbind(paste(selected_populations, collapse = ", "),nrow(my_Y_train), iteration_reliability))
    iteration_res_raw <- rbind(iteration_res_raw, iteration_pop_res)
    
    print(iteration_res_raw)
    
    if (iteration_pop_res$R2_testing <= min(iteration_res_raw$R2_testing)){ GAPIT_print_resplots(myGAPIT) }
  
  }
  
}

iteration_pop_sum <- iteration_res_raw %>% group_by (iteration) %>% summarise (rel_mean = mean(V2), rel_sd = sd(V2), R2_mean = mean(R2), R2_sd = sd(R2), R2_test_mean = mean(R2_testing), R2_test_sd = sd(R2_testing))


# 2. Validating by FID
FID_all <- unique(seedlings_all_library$NSWID)
runs = 5
iteration_FIDres_raw = NULL

for (iteration in iterations) {
  for (i in 1:runs){
    selected_FID <- sample(FID_all, 35)
    
    print(paste("Running pop:", paste(selected_FID, collapse = ", "), ". On run number", i, "out of", runs)) 
    
    my_Y_train <- Y[!(Y$NSWID %in% selected_FID),1:2]
    iteration_FID_res = NULL
    
    GAPIT_run_res <- GAPIT_function(iteration)
    
    myGAPIT <- GAPIT_run_res[[2]]
    iteration_reliability <- GAPIT_run_res[[1]]
    
    iteration_FID_res <- as.data.frame(cbind(paste(selected_FID, collapse = ", "),nrow(my_Y_train), iteration_reliability))
    iteration_FIDres_raw <- rbind(iteration_FIDres_raw, iteration_FID_res)
    
    print(iteration_FIDres_raw)
    
    if (iteration_FID_res$R2_testing <= min(iteration_FIDres_raw$R2_testing)){ GAPIT_print_resplots(myGAPIT) }
    
  }
  
}

iteration_FID_sum <- iteration_FIDres_raw %>% group_by (iteration) %>% summarise (rel_mean = mean(V2), rel_sd = sd(V2), R2_mean = mean(R2), R2_sd = sd(R2), R2_test_mean = mean(R2_testing), R2_test_sd = sd(R2_testing))

iteration_res_raw
iteration_pop_sum

iteration_FIDres_raw
iteration_FID_sum
