# Config values
setwd("~/ResearchCode/SNPsR")

gene_grouped_file <- "result/grouping-gene/10.GroupingComplete-custom.csv"
conclude_data <- read.csv(gene_grouped_file)

outfile_conclusion <- "result/skat/conclusion-p-value/conclusion.csv"

skat_result_path <- "result/skat/5-gene-test-new-0.0"
filez <- list.files(skat_result_path)

file_count <- 0
for (filename in filez) {
     p_value_file_path <- paste(skat_result_path, "/", filename, sep = "")
     p_value_data <- read.csv(p_value_file_path, header = TRUE)
     conclude_data <- cbind(conclude_data, p_value_data[, 6])
     
     file_count <- file_count + 1
}

colnames(conclude_data) <- c("Group No.",
                             "Group Name",	
                             "Start No.",
                             "Stop No.",
                             "Start BP",	
                             "Stop BP",
                             "New Members",
                             "Median",
                             paste("P", 1:file_count, sep = ""))

# Write all data to CSV
write.csv(conclude_data, file = outfile_conclusion, row.names = FALSE)

