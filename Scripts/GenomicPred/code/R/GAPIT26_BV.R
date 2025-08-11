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

iterations<- c(26)

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
  
  iteration_reliability <- rbind(iteration_reliability, cbind(iteration, mean(prediction_gt$Reliability), R2))
  ggsave (p2, file = "PredictionvPEV.jpg", height = 6.5, width = 6.5)
}


prediction_gt$Prediction[prediction_gt$Prediction < 0] <- 0

iteration_reliability

# Running blups without genotyped individual
library(tidyverse)
library(breedR)

# Joining seedling ID and library by IID
seedlings_all <- read.csv ("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Mquin_samples_pheno.csv")
seedlings_all_library <- merge(IID_height_sample_lib, seedlings_all, all.x=T, all.y=T)

genotyped_indivs <- myG[1,]
seedlings_all <- subset(seedlings_all_library, !(seedlings_all_library$LIBRARY %in% genotyped_indivs))
seedlings_all$COI <- seedlings_all$COI^0.25
seedlings_all <- seedlings_all %>% filter(!is.na(FID))

# Converting Mum FID to numeric
#NSWID <- unique(seedlings_all$IID)
mum_num <- max(seedlings_all$IID) + as.numeric(seq(1, length(unique(seedlings_all$FID))))
mum_num_fid <- as.data.frame(cbind(mum_num, unique(seedlings_all$FID)))
colnames(mum_num_fid)[2] <- "FID"
dataset <- left_join(seedlings_all, mum_num_fid)

dataset <- dataset %>% select (IID, mum_num, COI) %>% filter (!is.na(COI)) %>% mutate (dad = 0)
colnames(dataset) <- c('self', 'mum', "phe_X", 'dad')
dataset <- dataset [, c('self', 'dad', 'mum', 'phe_X')]

dataset$mum <- as.numeric(dataset$mum)


res.animal <- remlf90(fixed = phe_X ~ 1,
                      genetic = list(model = 'add_animal',
                                     pedigree = dataset[, 1:3],
                                     id = 'mum'),
                      data = dataset)
summary(res.animal)

## Extracting Predicted Breeding Values
## Predicted Breeding Values
# for the full pedigree first, 
PBV.full <- ranef(res.animal)$genetic #ranef extracts the conditional modes of the random effects from a fitted model object. For linear mixed models the conditional modes of the random effects are also the conditional means

# and for the observed individuals by matrix multiplication with the incidence matrix
PBV <- model.matrix(res.animal)$genetic %*% PBV.full

## Predicted genetic values vs. phenotype.
# Note: fitted = mu + PBV
qplot(fitted(res.animal), phe_X,
      data = dataset, colour = as.character(mum)) +
  geom_abline(intercept = 0,
              slope = 1,
              col = 'gray') +
  theme(legend.position='none')

fitted.res.animal.df <- data.frame(fitted(res.animal)) 

pred_gt_dataset <- cbind(dataset, fitted.res.animal.df)

# Plotting family BV
unique_mum <- unique(pred_gt_dataset$mum)
random_colours <- replicate(length(unique(pred_gt_dataset$mum)), paste0("#", paste(sample(0:9, 6, replace = TRUE), collapse = "")))
mum_colours <- data.frame(cbind(as.numeric(unique_mum), random_colours)) 
colnames(mum_colours)[1] <- "mum"
mum_colours$mum <- as.numeric(mum_colours$mum)

pred_gt_dataset_coloured <- left_join(pred_gt_dataset, mum_colours, copy = TRUE, by ="mum")

ggplot() + 
  geom_point(data = pred_gt_dataset_coloured, aes(x=mum, y=phe_X, color = random_colours), alpha = 0.3) +
  geom_point(data = pred_gt_dataset_coloured, aes(x=mum, y=fitted.res.animal., color = random_colours), shape=15, size = 2.5) +
  theme_bw() +
  scale_color_manual(values = pred_gt_dataset_coloured$random_colours, labels=pred_gt_dataset_coloured$mum) +
  theme(legend.position = "none")

# Generating GWAS phenotype dataset
parent_meta_subset <- parent_meta[,c(1,6)]
colnames(parent_meta_subset)[1] <- "FID"
pred_gt_dataset_subset <- unique(pred_gt_dataset[,c(3,5)])
mum_num_fid$FID <- as.character(mum_num_fid$FID)
parent_meta_subset$FID <- as.character(parent_meta_subset$FID)
BV_parent <- left_join(parent_meta_subset, mum_num_fid)
colnames(BV_parent)[3] <- "mum"
BV_parent$mum <- as.numeric(BV_parent$mum)
BV_parent <- left_join(BV_parent, pred_gt_dataset_subset, relationship = "many-to-many")
BV_parent <- BV_parent %>% filter(!is.na(fitted.res.animal.)) %>% select(sample_lib_NSW, fitted.res.animal.)


# Joining with blups

colnames(parent_meta_subset)[2] <- "Taxa"

prediction_gt_FID <- left_join(parent_meta_subset, prediction_gt)

colnames(BV_parent)[1] <- "Taxa"
BV_parent <- BV_parent %>% filter(!is.na(Taxa))
BV_valid <- left_join(BV_parent, prediction_gt_FID)
BV_valid_plot <- BV_valid %>% select (Taxa, fitted.res.animal., Prediction)
BV_valid_plot <- unique(BV_valid_plot)

ggplot(BV_valid_plot, aes(x=fitted.res.animal., y=Prediction)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (y = "Genomic Predicted COI^0.25", x = "Breeding-value COI^0.25") +
  theme_bw() +
  theme(legend.position = "right") +
  geom_text(mapping=aes(x=1.9, y= 0.1, label = "R2 = 0.24"), colour = "grey") +
  scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))

lm <- lm(Prediction ~ fitted.res.animal., data = subset(BV_valid_plot))
R2 <- summary(lm)$adj.r.squared
