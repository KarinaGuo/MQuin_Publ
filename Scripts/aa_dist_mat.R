library("Biostrings"); library("seqinr"); library("RColorBrewer"); library(gplots)
setwd("~/RBGSyd_Technical Officer/MQuin")

g97_fa_phase1 <- read.fasta(file = "snpEff/all_indivs/mq.phase1_edit_reverse_g97.fa")
g97_fa_phase2 <- read.fasta(file = "snpEff/all_indivs/mq.phase2_edit_reverse_g97.fa")

g93_fa_phase1 <- read.fasta(file = "snpEff/all_indivs/mq.phase1_edit_reverse_g93.fa")
g93_fa_phase2 <- read.fasta(file = "snpEff/all_indivs/mq.phase2_edit_reverse_g93.fa")

g4997_fa_phase1 <- read.fasta(file = "snpEff/all_indivs/mq.phase1_edit_reverse_g4997.fa")
g4997_fa_phase2 <- read.fasta(file = "snpEff/all_indivs/mq.phase2_edit_reverse_g4997.fa")

g4997_fa_phase1 <- lapply(g4997_fa_phase1, function(x) {comp(x)})
g4997_fa_phase2 <- lapply(g4997_fa_phase2, function(x) {comp(x)})

g6020_fa_phase1 <- read.fasta(file = "snpEff/all_indivs/mq.phase1_edit_reverse_g6020.fa")
g6020_fa_phase2 <- read.fasta(file = "snpEff/all_indivs/mq.phase2_edit_reverse_g6020.fa")


unique_ser=NULL
for (gene in c("g97", "g93", "g4997", "g6020")){
  print(paste("Started running", gene))
  
  fasta_phase1_varname <- paste0(gene,"_fa_phase1")
  fasta_phase2_varname <- paste0(gene,"_fa_phase2")
  
  fasta_phase1 <- get(fasta_phase1_varname)
  fasta_phase2 <- get(fasta_phase2_varname)
  
  aa_phase1 <- lapply((fasta_phase1), function(x) { seqinr::translate(x) })
  aa_phase2 <- lapply((fasta_phase2), function(x) { seqinr::translate(x) })
  
  unique_ser <- data.frame(rbind(unique_ser, cbind(gene=gene, unique_phase2 = sum(!duplicated(aa_phase2)))))
  
  aa_phases <- append(aa_phase1, aa_phase2)
  
  aa_al <- as.alignment(nb = length(aa_phases), seq = aa_phases, nam = names(aa_phases))
  
  aa_dist <- dist.alignment(aa_al)
  aa_dist_mat <- as.matrix(aa_dist)
  
  hm_gene <- paste0("heatmap_", gene)
  
  pdf(file = paste0("snpEff/all_indivs/heatmap_allindivs",gene,".pdf"))
  hm <- heatmap.2(
    as.matrix(aa_dist_mat),  # Distance matrix
    #margins = c(11, 11),  # Adjust margins
    trace = "none",  # No trace lines
    col = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),  # Heatmap colors
    key = T,  # Show color key
    density.info = "none",
    main=gene)
  dev.off()
  
  print(paste("Finished running", gene))
  
}

write.table(unique_ser, file = "snpEff/all_indivs/unique_ser.txt", row.names=FALSE, quote=FALSE)


# Heatmap of unique

for (gene in c("g97", "g93", "g4997", "g6020")){
  print(paste("Started running", gene))
  
  fasta_phase1_varname <- paste0(gene,"_fa_phase1")
  fasta_phase2_varname <- paste0(gene,"_fa_phase2")
  
  fasta_phase1 <- get(fasta_phase1_varname)
  fasta_phase2 <- get(fasta_phase2_varname)
  
  names(fasta_phase1) <- paste0(names(fasta_phase1),"_phase1")
  names(fasta_phase2) <- paste0(names(fasta_phase2),"_phase2")
  
  aa_phase1 <- lapply((fasta_phase1), function(x) { seqinr::translate(x) })
  # Unique indivs
  unique_aa_phase1_bin <- (duplicated(aa_phase1))
  unique_aa_phase1 <- names(aa_phase1)[!unique_aa_phase1_bin]
  
  aa_phase2 <- lapply((fasta_phase2), function(x) { seqinr::translate(x) })
  # Unique indivs
  unique_aa_phase2_bin <- (duplicated(aa_phase2))
  unique_aa_phase2 <- names(aa_phase2)[!unique_aa_phase2_bin]
  
  aa_phases <- append(aa_phase1, aa_phase2)
  unqiue_aa_phases <- append(unique_aa_phase1, unique_aa_phase2)
  
  # Recreate heatmap
  unique_aa_phase <- aa_phases[!is.na(match(names(aa_phases), unqiue_aa_phases))] 
  unique_aa_al <- as.alignment(nb = length(unique_aa_phase), seq = unique_aa_phase, nam = names(unique_aa_phase))
  
  unique_aa_dist <- dist.alignment(unique_aa_al)
  unique_aa_distmat <- as.matrix(unique_aa_dist)

pdf(file = paste0("snpEff/all_indivs/heatmap_allindivs_unique",gene,".pdf"))
heatmap.2(
  as.matrix(unique_aa_distmat),  # Distance matrix
  margins = c(11, 11),  # Adjust margins
  trace = "none",  # No trace lines
  col = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),  # Heatmap colors
  key = F,  # Show color key
  density.info = "none",
  main=gene)
dev.off()

print(paste("Finished running", gene))

}