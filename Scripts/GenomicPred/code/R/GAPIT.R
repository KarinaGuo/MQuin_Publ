# Following demo in https://zzlab.net/GAPIT/gapit_tutorial_script.txt
library(GAPIT)


#Import files
myY <- mq_phenotype_height_COI
myY$COI <- (myY$COI)^0.25
#myKI <- read.table("KSN.txt", head = FALSE)
Y=myY[,1:2]  
my_Y_train <- Y[-c(1:50),]
my_Y_valid <- Y[c(1:50),]

setwd("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_NLR/exact/unthinned")
myG <- read.csv("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/Mq_filt_cat_NLR_exact_unthin_sort.hapmap.hmp.txt", sep = "\t", header = FALSE)
myGAPIT <- GAPIT(Y=my_Y_train, G=myG, PCA.total=3, model=c("gBLUP"))

prediction=myGAPIT$Pred
unique(prediction$RefInf)

library(tidyverse)
colnames(myY)[1] <- "Taxa"
prediction_gt <- left_join(prediction, myY)

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

#ggsave (p1, file = "PredictionvGT.jpg")

#jpeg("Residuals.jpeg")
plot(prediction_gt_filtNA$Prediction, PredvCOI_resid, ylab="Residuals", xlab="Predicted COI"); abline(0, 0) 
#dev.off()


# From GAPIT.Association.Optimum.gBLUP.COI, genetic variance is 0.36
ga = 0.36
prediction_gt$Reliability = prediction_gt$PEV/ga

mean(prediction_gt$Reliability)

p2 <- ggplot(prediction_gt, aes(y=Reliability, x=Prediction)) +
  geom_point(aes(colour = RefInf)) +
  geom_smooth (method = "lm", se = F, linewidth = 0.5) +
  labs (x = "Genomic Predicted COI^0.25", y = "Reliability") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_discrete(name = "Prediction Category", labels = c("Training", "Testing"))
p2

#ggsave (p2, file = "PredictionvPEV.jpg")



#### For GAPIT_6&7

# Replacing problematic Positions with those of relevant clusters. SNPs are saved in all_table_sigSNPs.txt, which the positions were gathered from iter_2 gwas_2 and the values calculated from /iter_4/data/catvcf/Mq_filt_cat.vcf.gz
# Relevant SNPs of the clusters are found in /home/karina/mqgwas/iter_5/gwas_01_impute/Output/Linear_Mixed_Model/. E.g, sig_SNPs_6_mod_sub_clustered_01POS_beagle_heterozygosity.part6_Mq_filt_cat.beagle.csv
all_SNPs <- read.table("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_7/all_table_sigSNPs.txt", header = TRUE)
all_SNPs$POS <- as.numeric(all_SNPs$POS)

# prob_SNP_POS_1 is the ideal ones, but we will ignore everything that was not used in the previous clustering and use POS_2
prob_SNP_POS_1 <- c(10567634, 10569230, 19224573, 19224562, 19224521, 19561425, 19291828, 19352903, 19353548, 19336837, 19354859, 19358240, 19348984, 19349112, 19349526, 19347499, 19347539, 19335439, 19335449, 19348177, 19348203, 19347642, 19349286, 19349317, 19339025, 19349473, 19519429) # selected manually by looking at values
prob_SNP_POS_2 <- c(19224573, 19224562, 19224521, 19561425, 19291828, 19352903, 19353548, 19336837, 19354859, 19358240, 19348984, 19349112, 19349526, 19347499, 19347539, 19335439, 19335449, 19348177, 19348203, 19347642, 19349286, 19349317, 19339025, 19349473, 19519429) # selected manually by looking at values
  
# Getting clusters from genotype_as_phenotype.R in iter_5
cluster_groups <- POS_01_grouped
colnames(cluster_groups)[1] <- "POS"

ggplot() + 
  geom_jitter(aes(y=POS_01_grouped_cols$POS, x=0, col = POS_01_grouped_cols$random_colors), width = 0.01, alpha=0.8) +
  geom_jitter(aes(y=all_SNPs$POS, x=0.05), width = 0.01, alpha=0.8) +
  geom_jitter(aes(y=prob_SNP_POS_2, x=0.05, col = "red"), width = 0.01, alpha=0.8) +
  #scale_y_continuous(breaks = pretty(POS_01_grouped_cols$POS, n = 30)) +
  ylim(19000000, 19600000) +
  xlim(-0.1,0.15) +
  theme_minimal() +
  theme(legend.position = "none") 

# Finding which problematic SNP is closest to which cluster 
cluster_mean <- cluster_groups %>% group_by (clust) %>% summarise (mean = mean(POS))

prob_POS_clust <- data.frame()
for (POS in prob_SNP_POS_2) {
  POS_dist <- as.numeric((POS-cluster_mean$mean)^2)
  clust_ind <- match (min(POS_dist), POS_dist)
  prob_POS_clust <- rbind(prob_POS_clust, cbind(POS, as.integer(clust_ind)))
}

paste("Getting results for", sort(unique(prob_POS_clust$V2)))
colnames(prob_POS_clust)[2] <- "clust"

# Reading in sigSNP results per cluster

SNP_repl_clusters <- data.frame()

for (repl_clust in c(unique(prob_POS_clust$clust))){
   fp <- paste0("~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_7/sig_SNPs_",repl_clust,"_mod_sub_clustered_01POS_beagle_heterozygosity.part",repl_clust,"_Mq_filt_cat.beagle.csv")
   SNP_repl_clust <- read.csv(fp)
   SNP_repl_clust <- SNP_repl_clust[,c(1,3)]
   SNP_repl_clust$clust <- repl_clust
   SNP_repl_clusters <- rbind(SNP_repl_clusters, SNP_repl_clust)
}
colnames(SNP_repl_clusters) <- c("CHROM", "POS", "clust")


all_SNPs_rmprob <- all_SNPs[!(all_SNPs$POS %in% prob_SNP_POS_2),]
all_SNPs_rmprob <- all_SNPs_rmprob[, c("CHROM", "POS")]
SNP_repl_clusters <- SNP_repl_clusters[, c("CHROM", "POS")]

all_SNPs_probrepl <- rbind(all_SNPs_rmprob, SNP_repl_clusters)
all_SNPs_probrepl <- unique(all_SNPs_probrepl)
write.table(all_SNPs_probrepl, file="~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_7/GAPIT_7_SNPs_replaced.txt", row.names=FALSE, quote = FALSE, col.names = FALSE)
write.table(all_SNPs_rmprob, file="~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/data/GAPIT_6/GAPIT_6_SNPs_rm.txt", row.names=FALSE, quote = FALSE, col.names = FALSE)

# GAPIT 25 Clustering populations, removing outlying populations as the testing pops
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

my_Y_train <- subset(Y, !(Y$LIBRARY %in% samples_testing$LIBRARY))
# nrow



### Plotting GAPIT quality v iteration (1 - 14)
GAPIT_res <- read.csv(file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/GAPIT_results_plot_231123.csv")


p1 <- ggplot() + 
  geom_point(data = GAPIT_res, mapping = aes (x=log(SNP.count), y = R2, colour = Type)) +
  geom_smooth(data = GAPIT_res, mapping = aes (x=log(SNP.count), y = R2), method = 'lm', se = FALSE, linewidth = 0.7, colour = "black") +
  ggrepel::geom_text_repel(data = GAPIT_res, mapping = aes (x=log(SNP.count), y = R2, colour = Type, label = Iteration),  min.segment.length = 0.005) +
  theme_bw()

ggsave (p1, file = "~/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/GAPIT_results_plot_231123.jpg", width = 2000, height = 1500, limitsize = FALSE, units = "px")

