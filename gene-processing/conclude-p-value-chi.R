# Config values
setwd("~/ResearchCode/SNPsR")

snpsName_path <- "input/gen-snps/NameONLY13479.csv"
conclude_data <- read.csv(snpsName_path, header = TRUE)

outfile_conclusion <- "result/skat/conclusion-p-value/conclusion-chi-square.csv"

chisquare_result_path <- "result/chi-squared/0.0"
filez <- list.files(chisquare_result_path)

file_count <- 0
for (filename in filez) {
     p_value_file_path <- paste(chisquare_result_path, "/", filename, sep = "")
     p_value_data <- read.csv(p_value_file_path, header = TRUE)
     conclude_data <- cbind(conclude_data, p_value_data[, 3])
     
     file_count <- file_count + 1
}

colnames(conclude_data) <- c("SNPs Name",
                             paste("P", 1:file_count, sep = ""))

# Write all data to CSV
write.csv(conclude_data, file = outfile_conclusion, row.names = FALSE)

