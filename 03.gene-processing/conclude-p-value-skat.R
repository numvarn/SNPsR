# Config values
setwd("~/ResearchCode/SNPsR")

gene_grouped_file <- "result/grouping-gene/10.GroupingComplete.csv"
conclude_data <- read.csv(gene_grouped_file)

outfile_conclusion <- "result/skat/conclusion-p-value/conclusion-skat-3000replicated.csv"

# Floder that stroe skat results
skat_result_path <- "result/skat/0.0/linear-weighted"
filez <- list.files(skat_result_path)

filename_list <- c()

file_count <- 0
for (filename in filez) {
     p_value_file_path <- paste(skat_result_path, "/", filename, sep = "")
     p_value_data <- read.csv(p_value_file_path, header = TRUE)
     conclude_data <- cbind(conclude_data, p_value_data[, 6])
     
     cat(sprintf("Processing filename : %s\n", filename))
     filename_list <- c(filename_list, filename)
}

colnames(conclude_data) <- c("Group No.",
                             "Group Name",	
                             "Start No.",
                             "Stop No.",
                             "Start BP",	
                             "Stop BP",
                             "New Members",
                             "Median",
                             paste("P-", filename_list, sep = ""))

# Write all data to CSV
write.csv(conclude_data, file = outfile_conclusion, row.names = FALSE)