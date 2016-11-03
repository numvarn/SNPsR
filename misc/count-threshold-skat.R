# Auther Phisan Sookkhee

setwd("~/ResearchCode/SNPsR")

conclude_skat_file <- "result/skat/conclusion-p-value/conclusion-skat-3000replicated.csv"
conclude_out_file <- "result/skat/conclusion-p-value/conclusion-skat3000threshold.csv"

cat(sprintf("Reding file %s\n", conclude_skat_file))
conclude_skat_data <- read.csv(conclude_skat_file, header = TRUE)

col_number <- ncol(conclude_skat_data)
row_number <- nrow(conclude_skat_data)

threshold_bonf <- -log10(0.05 / row_number)
threshold_bonf_norm <- -log10(0.05)

# Create matrix for store result
conclude_skat_mx <- matrix(NA, row_number+3, col_number)
for (i in 1:col_number) {
     if (i <= 8) {
          for (j in 1:row_number) {
               conclude_skat_mx[j, i] <- as.character(conclude_skat_data[j, i])
          }  
          if (i == 8) {
               conclude_skat_mx[row_number+1, i] <- as.character(paste(">-log10(0.05 / ", row_number, ")"))
               conclude_skat_mx[row_number+2, i] <- as.character(">-log10(0.05)")
               conclude_skat_mx[row_number+3, i] <- as.character("status")
          }
     } else {
          cat(sprintf("Processing SNPs #%s\n", i-8))
          threshold_bonf_count <- 0
          threshold_bonf_norm_count <- 0
          status <- 0
          
          for (j in 1:row_number) {
               conclude_skat_mx[j, i] <- conclude_skat_data[j, i]
               
               if (conclude_skat_data[j, i] > threshold_bonf_norm) {
                    threshold_bonf_norm_count <- threshold_bonf_norm_count + 1
               }
               
               if (conclude_skat_data[j, i] > threshold_bonf) {
                    threshold_bonf_count <- threshold_bonf_count + 1
                    status <- 1
               }
          }
          conclude_skat_mx[row_number+1, i] <- threshold_bonf_count
          conclude_skat_mx[row_number+2, i] <- threshold_bonf_norm_count
          conclude_skat_mx[row_number+3, i] <- status
     }
}

cat(sprintf("Writing result to %s\n", conclude_out_file))
write.csv(conclude_skat_mx, file = conclude_out_file, row.names = FALSE)
