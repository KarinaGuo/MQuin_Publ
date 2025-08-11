# run Calculating metrics.R which is used for the metadata

seedlings_all <- read.csv ("~/RBGSyd_Technical Officer/MQuin/Processing Meta/Mquin_samples_pheno.csv") # Meta data for seedlings including height and COI of all seedlings, including those that didn't get genotypes
colnames(seedlings_all)[which(names(seedlings_all)=="FID")] <- "NSWID"
library(tidyverse)
library(breedR)
library(lme4)


### For raw COI
# Converting Mum FID to numeric
NSWID <- unique(seedlings_all$NSWID)
mum_num <- max(seedlings_all$IID) + as.numeric(seq(1, length(NSWID)))
mum_num_fid <- as.data.frame(cbind(mum_num, NSWID))
dataset <- merge(seedlings_all, mum_num_fid, by.x ='NSWID', by.y='NSWID')

dad_num <- seq(from=max(mum_num),by=1,to=as.numeric(max(mum_num)+nrow(dataset)-1))
dataset <- dataset %>% dplyr::select (IID, mum_num, COI) %>% mutate (dad = dad_num) %>% filter (!is.na(COI)) 
colnames(dataset) <- c('self', 'mum', "phe_X", 'dad')
dataset <- dataset [, c('self', 'dad', 'mum', 'phe_X')]


dataset$mum <- as.numeric(dataset$mum)

#dataset$phe_X <- dataset$phe_X^0.25
dataset$phe_X <- dataset$phe_X
dataset <- dataset %>% group_by (mum) %>% mutate(n=n()) %>% ungroup() %>% filter(n>3) %>% select(-c(n)) %>% data.frame() # 1020 indivs vs orig 1098


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
fitted.res.animal.df$index <- seq(1, nrow(fitted.res.animal.df), 1)

dataset.index <- dataset
dataset.index$index <- seq(1, nrow(dataset), 1)

pred_gt_dataset <- left_join(dataset.index, fitted.res.animal.df)

# Plotting family BV
unique_mum <- unique(pred_gt_dataset$mum)
random_colours <- replicate(length(unique(pred_gt_dataset$mum)), {
  rgb <- sample(50:255, 3, replace = TRUE) # Higher range for brighter colors
  paste0("#", sprintf("%02X%02X%02X", rgb[1], rgb[2], rgb[3]))
})
mum_colours <- data.frame(cbind(as.numeric(unique_mum), random_colours)) 
colnames(mum_colours)[1] <- "mum"
mum_colours$mum <- as.numeric(mum_colours$mum)

pred_gt_dataset_coloured <- left_join(pred_gt_dataset, mum_colours, copy = TRUE, by ="mum")
ordered_pred_gt_dataset_coloured <- pred_gt_dataset_coloured %>%
  arrange(fitted.res.animal.)
ordered_pred_gt_dataset_coloured$mum <- factor(
  ordered_pred_gt_dataset_coloured$mum,
  levels = unique(ordered_pred_gt_dataset_coloured$mum)
)

## Removing individual with only one sample
ordered_coloured_rmlow <- ordered_pred_gt_dataset_coloured %>% group_by(mum) %>% mutate (n=n()) %>% ungroup() %>% filter(n>2)

plotly::ggplotly(ggplot() + 
  geom_point(data = ordered_pred_gt_dataset_coloured, aes(x=mum, y=phe_X, color = random_colours), alpha = 0.3) +
  geom_boxplot(data = ordered_pred_gt_dataset_coloured, aes(x=mum, y=phe_X, color = random_colours), alpha = 0.1) +
  geom_point(data = ordered_pred_gt_dataset_coloured, aes(x=mum, y=fitted.res.animal.),color = "black", shape=15, size = 1) +
  theme_bw() +
  #scale_color_manual(values = pred_gt_dataset_coloured$random_colours, labels=pred_gt_dataset_coloured$mum) +
  #theme(legend.position = "none",  axis.text.x=element_blank(),  axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 20)) +
  theme(legend.position = "none") +
  labs(y="COI", x = "Maternal ID (Sorted)"))

# Generating GWAS phenotype dataset
parent_meta_subset <- parent_meta[,c("NSWID","sample_lib_NSW")]
colnames(parent_meta_subset)[1] <- "NSWID"
pred_gt_dataset_subset <- unique(pred_gt_dataset[,c(3,6)])
BV_parent <- left_join(parent_meta_subset, mum_num_fid)
colnames(BV_parent)[3] <- "mum"
BV_parent$mum <- as.numeric(BV_parent$mum)
BV_parent <- merge(BV_parent, pred_gt_dataset_subset)

write.csv(BV_parent, file ="~/RBGSyd_Technical Officer/MQuin/Parent/Iteration 1/data/metabreedingcoeff_parent_animalmdlBV_2.csv",row.names=FALSE)


BV_parent_sort <- BV_parent[order(BV_parent$fitted.res.animal.),] %>% filter(!is.na(fitted.res.animal.) & !is.na(sample_lib_NSW))
BV_parent_sort_uniq <- BV_parent_sort[!duplicated(BV_parent_sort[,c("mum", "fitted.res.animal.")]),]
BV_parent_sort_uniq[c(1:10, (nrow(BV_parent_sort_uniq)-10):(nrow(BV_parent_sort_uniq))),]

write.table(BV_parent_sort_uniq[c(1:10, (nrow(BV_parent_sort_uniq)-10):(nrow(BV_parent_sort_uniq))),c("sample_lib_NSW")], file = "~/RBGSyd_Technical Officer/MQuin/Parent/SMCPP/Run_11/indvs.txt", row.names = F, col.names = F, quote = F)

################################################################## For raw RGR
# Calculating RGR
sample_list_mquin_PBI$HT1 <- as.numeric(sample_list_mquin_PBI$HT1)
sample_list_mquin_PBI$HT2 <- as.numeric(sample_list_mquin_PBI$HT2)
#Calculating RGR
sample_list_mquin_PBI$RGR <- (log(sample_list_mquin_PBI$HT2) - log(sample_list_mquin_PBI$HT1))/25 #Date between 23/08/23 and 16/09/23)
sample_list_mquin_PBI$RGR[sample_list_mquin_PBI$HT2 == 0] <- 0

seedlings_all <- left_join(seedlings_all, sample_list_mquin_PBI)

# Converting Mum FID to numeric
NSWID <- unique(seedlings_all$NSWID)
mum_num <- max(seedlings_all$IID) + as.numeric(seq(1, length(NSWID)))
mum_num_fid <- as.data.frame(cbind(mum_num, NSWID))
dataset <- left_join(seedlings_all, mum_num_fid)

dataset <- dataset %>% select (IID, mum_num, RGR) %>% filter (!is.na(RGR)) %>% mutate (dad = 0)
colnames(dataset) <- c('self', 'mum', "phe_X", 'dad')
dataset <- dataset [, c('self', 'dad', 'mum', 'phe_X')]

dataset$mum <- as.numeric(dataset$mum)


res.animal_RGR <- remlf90(fixed = phe_X ~ 1,
                      genetic = list(model = 'add_animal',
                                     pedigree = dataset[, 1:3],
                                     id = 'mum'),
                      data = dataset)
summary(res.animal_RGR)

## Extracting Predicted Breeding Values
## Predicted Breeding Values
# for the full pedigree first, 
PBV.full <- ranef(res.animal_RGR)$genetic #ranef extracts the conditional modes of the random effects from a fitted model object. For linear mixed models the conditional modes of the random effects are also the conditional means

# and for the observed individuals by matrix multiplication with the incidence matrix
PBV <- model.matrix(res.animal_RGR)$genetic %*% PBV.full

## Predicted genetic values vs. phenotype.
# Note: fitted = mu + PBV
qplot(fitted(res.animal_RGR), phe_X,
      data = dataset) +
  geom_abline(intercept = 0,
              slope = 1,
              col = 'gray')

fitted.res.animal.df <- data.frame(fitted(res.animal_RGR)) 
fitted.res.animal.df$index <- seq(1, nrow(fitted.res.animal.df), 1)

dataset.index <- dataset
dataset.index$index <- seq(1, nrow(dataset), 1)

pred_gt_dataset <- left_join(dataset.index, fitted.res.animal.df)

# Plotting family BV
unique_mum <- unique(pred_gt_dataset$mum)
random_colours <- replicate(length(unique(pred_gt_dataset$mum)), paste0("#", paste(sample(0:9, 6, replace = TRUE), collapse = "")))
mum_colours <- data.frame(cbind(as.numeric(unique_mum), random_colours)) 
colnames(mum_colours)[1] <- "mum"
mum_colours$mum <- as.numeric(mum_colours$mum)

pred_gt_dataset_coloured <- left_join(pred_gt_dataset, mum_colours, copy = TRUE, by ="mum")

ggplot() + 
  geom_point(data = pred_gt_dataset_coloured, aes(x=mum, y=phe_X, color = random_colours), alpha = 0.3) +
  geom_point(data = pred_gt_dataset_coloured, aes(x=mum, y=fitted.res.animal_RGR., color = random_colours), shape=15, size = 2.5) +
  theme_bw() +
  scale_color_manual(values = pred_gt_dataset_coloured$random_colours, labels=pred_gt_dataset_coloured$mum) +
  theme (legend.position = "none")

# Generating GWAS phenotype dataset
parent_meta_subset <- parent_meta[,c(5,6)]
colnames(parent_meta_subset)[1] <- "NSWID"
pred_gt_dataset_subset <- unique(pred_gt_dataset[,c(3,6)])
BV_parent_RGR <- left_join(parent_meta_subset, mum_num_fid)
colnames(BV_parent_RGR)[3] <- "mum"
BV_parent_RGR$mum <- as.numeric(BV_parent_RGR$mum)
BV_parent_RGR <- left_join(BV_parent_RGR, pred_gt_dataset_subset, relationship = "many-to-many")
BV_parent_RGR <- BV_parent_RGR %>% filter(!is.na(fitted.res.animal_RGR.)) %>% select(sample_lib_NSW, fitted.res.animal_RGR.)

BV_parent_both <- left_join(BV_parent, BV_parent_RGR)
# ggplot(BV_parent_both, aes(x=fitted.res.animal., y=fitted.res.animal_RGR.)) + geom_point() + geom_smooth(method = "lm", se=F) + theme_bw()

#write.csv(BV_parent_both, file ="~/RBGSyd_Technical Officer/MQuin/Parent/data/meta/breedingcoeff_parent_animalmdlBV_RGR.csv",row.names=FALSE)