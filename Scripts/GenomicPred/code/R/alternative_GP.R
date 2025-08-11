library(BGLR)
library(BGData)

################## Attemppt 1
write.table(demo(BRR), file = "~/Delete/BGLR_Demo.txt")
rm(list=ls())
setwd(tempdir())
library(BGLR)
data(wheat) #Loads the wheat dataset Matrix Y contains the average grain yield, column 1: Grain yield for environment 1 and so on. The matrix A contains additive relationship computed from the pedigree and matrix X contains the markers information. 
# X is a matrix (599 x 1279) with DArT genotypes; data are from pure lines and genotypes were coded as 0/1 denoting the absence/presence of the DArT. Markers with a minor allele frequency lower than 0.05 were removed, and missing genotypes were imputed with samples from the marginal distribution of marker genotypes, that is, x_{ij}=Bernoulli(\hat p_j), where \hat p_j is the estimated allele frequency computed from the non-missing genotypes. 

y=wheat.Y[,1]
### Creates a testing set with 100 observations
whichNa<-sample(1:length(y),size=100,replace=FALSE)
yNa<-y
yNa[whichNa]<-NA
### Runs the Gibbs sampler
fm<-BLR(y=yNa,XL=wheat.X,GF=list(ID=1:nrow(wheat.A),A=wheat.A),
        prior=list(varE=list(df=3,S=0.25),
                   varU=list(df=3,S=0.63),
                   lambda=list(shape=0.52,rate=1e-4,
                               type='random',value=30)),
        nIter=500,burnIn=100,thin=1)

#
# y is the phenotype
NSWID <- unique(seedlings_all$NSWID)
mum_num <- max(seedlings_all$IID) + as.numeric(seq(1, length(NSWID)))
mum_num_fid <- as.data.frame(cbind(mum_num, NSWID))
dataset <- left_join(seedlings_all, mum_num_fid)

dataset <- dataset %>% select (IID, mum_num, COI) %>% mutate (dad = 0)
colnames(dataset) <- c('self', 'mum', "phe_X", 'dad' )
dataset <- dataset [, c('self', 'dad', 'mum', 'phe_X')]

dataset$mum <- as.numeric(dataset$mum)

seedling.Y <- dataset[,c(1,4)]

# Calculating A
library(optiSel)

pedigree <- dataset[,-4]

seedling.A <- makeA(pedigree, AFounder=NULL)

# Converting to genotype x (presence/absence marker matrix)


################# Attempt 2
#1. Flat Prior (FIXED)

library(BGLR)

pheno <- seedling.Y
attach(pheno)

X <- read.table(file="C:/Users/swirl/OneDrive/Documents/RBGSyd_Technical Officer/MQuin/Seedling GWAS/Filtering/Iteration 2/data/evalmetrics/evalmetric_Mq_GT_sigSNP.out.GT.FORMAT")

fm=BGLR(y=

data(mice); X=scale(mice.X[,1:3000]); pheno=mice.pheno


fm=BGLR(y=Obesity.BMI,ETA=list( list(~GENDER+CoatColour+CageDensity,model='FIXED')), nIter=6000,burnIn=1000)
fm2=lm(Obesity.BMI~GENDER+CoatColour+CageDensity)
plot(cbind(c(fm$mu,fm$ETA[[1]]$b),coef(fm2))); abline(a=0,b=1)
