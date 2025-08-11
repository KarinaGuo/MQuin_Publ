GAPIT_res <- read.csv(file = "C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 6/GAPIT_results_plot_240603.csv")

library(tidyverse)
library(ggrepel)
p1 <- ggplot(GAPIT_res, aes(x=log(SNP.count), y=R2)) +
  geom_point(aes(colour=Type)) +
  theme_bw() +
  geom_text_repel(aes(label=Iteration), min.segment.length=0) +
  labs(x="Logged number of SNPs", y="R2")
  
p1

p2 <- ggplot(GAPIT_res, aes(x=log(SNP.count), y=R2)) +
  stat_smooth(colour='grey', se=F, size=1, method='lm') +
  geom_point(aes(colour=Type)) +
  facet_wrap(vars(Type), nrow=2, scales = "free_x") +
  theme_bw() +
  geom_text_repel (aes(label=Iteration), size=4, min.segment.length=0) +
  labs(x="Logged number of SNPs", y="R2")
p2

GAPIT_res_markers <- GAPIT_res %>% filter(Type=="Markers") 

p3 <- ggplot(GAPIT_res_markers, aes(x=SNP.count, y=Crossval_R2)) +
  stat_smooth(colour='grey', se=F, size=1, method='lm') +
  geom_point() +
  theme_bw() +
  geom_text_repel(label=GAPIT_res_markers$Iteration, min.segment.length=0) +
  labs(x="Number of SNPs", y="R2 of only cross-validation individuals")
  
p3
