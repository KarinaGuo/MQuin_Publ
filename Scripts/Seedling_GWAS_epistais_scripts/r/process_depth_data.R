#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

infile <- paste0(args[1], "/data/evalmetric_Mq_depth_phenosub_file.out.idepth")
outfile <- paste0(args[1], "/data/meta/samples_gt_depth_thresh.txt")
qd <- read.table(infile, sep="\t",header=TRUE)

calculate_weighted_mean <- function(qd) {
  grouped <- split(qd, qd$INDV)

  weighted_means <- sapply(grouped, function(group) {
    weighted.mean(group$MEAN_DEPTH, group$N_SITES)
  })

  result <- data.frame(INDV = names(weighted_means), Weighted_Mean = weighted_means)
  return(result)
}

weighted_means <- calculate_weighted_mean(data)

i_gt_thresh <- which(qd$Weighted_Mean >= 3)
write.table(qd$INDV[i_gt_thresh], outfile,col.names=FALSE, row.names=FALSE, quote=FALSE), quote=FALSE)