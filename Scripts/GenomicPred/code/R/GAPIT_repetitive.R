# GAPIT Iterative validation
library(GAPIT)
library(tidyverse)

repitions = 5
test_percentage = 1/10 #(52)

myY <- mq_phenotype_height_COI %>% filter (!is.na(COI))
colnames(myY) [1] <- "Taxa"
Y=myY[,1:2] # traits
Y[,2] <- Y[,2]^0.25

n=nrow(Y) # Number of total
testing=round(n*test_percentage) # Training

storage=NULL

chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)

iterations <- c(23)

iteration_reliability <- data.frame()

for (repetition in iterations){
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

  ## Iteratively running results on random testing samples, number of rounds = repetitions

  for(rep in 1:repitions){
    #Set missing data
    sample.missing=sample(1:n,testing)
    Y0=Y[-sample.missing,]
    #Prediction
    remove(myGAPIT)
    myGAPIT <- tryCatch(expr = {(GAPIT(Y=Y0, G=myG, PCA.total=3, model=c("gBLUP")))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=Y0, G=myG, PCA.total=3, model=c("gBLUP"), Geno.View.output = FALSE))})
    
    prediction=myGAPIT$Pred
    
    #Separate reference (with phenotype) and inference (without phenotype)
    prediction.ref=prediction[prediction[,3]==1,]
    prediction.inf=prediction[prediction[,3]==2,]
    
    #Merge prediction with original Y
    YP.ref <- merge(myY, prediction.ref, by.x = "Taxa", by.y = "Taxa")
    YP.inf <- merge(myY, prediction.inf, by.x = "Taxa", by.y = "Taxa")
    
    #Calculate correlation and store them
    r.ref=(cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,10]) )))^2
    r.inf=(cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,10]) )))^2
    
    repetition_ref <- cbind(iteration, rep, r.ref)
    repetition_inf <- r.inf
    
    storage=as.data.frame(rbind(storage, cbind(repetition_ref,repetition_inf)))
  }#End of for (rep in 1:t)

colnames(storage)=c("Iteration", "Repetition", "Reference_R2","Inference_R2")
iteration_summ <- storage %>% group_by (Iteration) %>% summarise (mean_Reference_R2 = mean(Reference_R2), sd_Reference_R2 = sd(Reference_R2), mean_Inference_R2 = mean(Inference_R2), sd_Inference_R2 = sd(Inference_R2) )
