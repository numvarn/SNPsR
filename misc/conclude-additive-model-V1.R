# Code by Phisan Sookkhee
# Update 27 Dec 2016

library("gtools")

# Config values
setwd("~/ResearchCode/SNPsR")

snpsName_path <- "input/gen-snps/NameONLY13479.csv"
conclude_data <- read.csv(snpsName_path, header = TRUE)

outfile_conclusion <- "/Users/phisan/Dropbox/04.Research/NumvarnThesis/conclud_and_count/results/result-conclude.csv"

glm_result_path <- "/Users/phisan/Dropbox/04.Research/NumvarnThesis/conclud_and_count/data"
filez <- list.files(glm_result_path)
filez <- mixedsort(filez)

row_number <- nrow(conclude_data)
error_list <- c()

file_count <- 0
for (filename in filez) {
     p_value_file_path <- paste(glm_result_path, "/", filename, sep = "")
     p_value_data <- read.csv(p_value_file_path, header = TRUE)
     
     if(nrow(p_value_data) == row_number) {
          conclude_data <- cbind(conclude_data, p_value_data[, 3])
     } else {
          error_list <- c(error_list, filename)
     }
     
     file_count <- file_count + 1
     
     cat(sprintf("Procssing file no. #%s : %s\n", file_count, filename))
}

# Show error file 
if(length(error_list) > 0) {
     cat(sprintf("Error file is "))
     for(filename in error_list) {
          cat(sprintf("error : %s\n", filename))
     }
}

colnames(conclude_data) <- c("SNPs Name",
                             paste("P", 1:file_count, sep = ""))

# Write all data to CSV
write.csv(conclude_data, file = outfile_conclusion, row.names = FALSE)

