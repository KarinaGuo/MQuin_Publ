# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#library(GAPIT)
library(tidyverse)

#Import files
myY <- mq_phenotype_height_COI
myY$COI <- (myY$COI)^0.25
Y=myY[,1:2]  
my_Y_train <- Y[-c(1:150),]
# Replacing Chrom names
chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

GWAS_res <- read.csv(file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 4/gwas/gwas_2/Output/Linear_Mixed_Model/COI/COI_20231025_215412/best_p-values/p_wald_filenamecut_top0.01.csv")
genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_40_sigSNP.hapmap.hmp.txt")
myG <- read.csv(file = genetic_data, sep = "\t", header = FALSE)

# Removing which POS
myG_sub <- myG %>% dplyr::select(V3, V4); myG_sub <- myG_sub[-c(1),];  colnames(myG_sub) <- c("chr", "ps")
GWAS_res$ps <- as.character(GWAS_res$ps); GWAS_res$chr <- gsub("Mq", "MQ", GWAS_res$chr)
myG_score <- left_join(myG_sub, GWAS_res)

myG_ordered <- myG_score[order(myG_score$p_wald,decreasing = T),]
myG_to_rm <- myG_ordered[1:50,1:2]

# Removing POS
myG_to_rm[1,] <- names(myG_to_rm) #colnames(myG) <- names(myG_to_rm)
myg40 <- myG %>% filter(!(myG$V4 %in% myG_to_rm$ps))
nrow(myG) - nrow(myg40)

# Adding Top 50 pos epi hits
epistatic_hits <- read.table("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 4/Epistasis/whole_significant/plink_whole_full_noepi1_epistasis_plinked.epi.cc.summary.MRsigSNPfilt")
epistatic_hits_ord <- epistatic_hits[order(epistatic_hits$V6, decreasing = T),]; epistatic_hits_ord_top50 <- epistatic_hits_ord[50,]

## Top 50 
epistatic_hits_ord_split_data <- strsplit(as.character(epistatic_hits_ord$V2), ":")
epistatic_hits_ord_split <- as.data.frame(do.call(rbind, epistatic_hits_ord_split_data))
rm_MR <- epistatic_hits_ord_split %>% filter (!(V2 %in% myG$V4)) # removing epistatic hits that are already MR markers
rm_MR_top50 <- rm_MR[1:50,] # top 50

# Subset top 50
genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_42_sigSNP.hapmap.hmp.txt")
myG_42 <- read.csv(file = genetic_data, sep = "\t", header = FALSE)

myG_42_top50 <- myG_42 %>% filter(V4 %in% rm_MR_top50$V2)

# Binding
myg43 <- rbind(myg40,myG_42_top50 )
nrow(myg43)
#write.table(myg43, sep = "\t", row.names = FALSE, col.names = FALSE, file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_43_sigSNP.hapmap.hmp.txt")



## Running
iteration <- 43
iteration_reliability <- data.frame()

file_dir = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration)
wd = file_dir
genetic_data = paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_",iteration,"_sigSNP.hapmap.hmp.txt")

# Running gBLUP for GP
setwd(wd)
myG <- read.csv(file = genetic_data, sep = "\t", header = FALSE)
myG$V3[1] <- "chrom"
myG$V3 <- chromosome_mapping[myG$V3] 
if(!file.exists(file_dir)){
  print(paste("Directory", file_dir, "does not exist and is being created :)")) 
  dir.create(file_dir)
}


# Running gBLUP for GP
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

p1 <- ggplot(prediction_gt[prediction_gt$RefInf==c(1,2),], aes(x=COI, y=Prediction)) +
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

lm <- lm(Prediction ~ COI, data = subset(prediction_gt, RefInf==2))
R2_testing <- summary(lm)$adj.r.squared

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

iteration_reliability <- rbind(iteration_reliability, cbind(iteration, mean(prediction_gt$Reliability), R2, R2_testing))
ggsave (p2, file = "PredictionvPEV.jpg", height = 6.5, width = 6.5)

iteration_reliability

