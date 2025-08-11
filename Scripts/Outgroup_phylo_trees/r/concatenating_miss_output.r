directory <- "/home/karina/mqgwas/outgroups_iter_1/data/BUSCO/evalmetric_GQ20_BUSCO_missing/"
file_list <- list.files(directory)

# Creating an initial dataframe with the first file
i_first <- file_list[1]
filepath = paste0(directory,i_first)
file <- read.table(filepath, header = TRUE)
file$N_PROP <- file$N_MISS/file$N_DATA
file <- file [,c(1,6)]
colnames(file)[2] <- i_first

depth_concat = file

# Removing first occurrence as used above
file_list <- file_list[-1]

# Repeating process and left joining
for (i in file_list){
  filepath = paste0(directory,i)
  file_name <- i
  file <- read.table(filepath, header = TRUE)
  file$N_PROP <- file$N_MISS/file$N_DATA
  file <- file [,c(1,6)]
  colnames(file)[2] <- file_name
  depth_concat  <- merge(x = depth_concat, y = file, by = "INDV", all.x=TRUE)
}

out_file <- paste0(directory, "concatenated_results.csv")
write.csv (depth_concat, file = out_file)