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

# Create matrix for store result.
conclude_skat_mx <- matrix(0, col_number-7, 4)

for (i in 9:col_number) {
     index <- i - 8
     cat(sprintf("Processing Replicated #%s\n", index))
     threshold_bonf_count <- 0
     threshold_bonf_norm_count <- 0
     status <- 0
     
     for (j in 1:row_number) {
          if (conclude_skat_data[j, i] > threshold_bonf_norm) {
               threshold_bonf_norm_count <- threshold_bonf_norm_count + 1
          }
          
          if (conclude_skat_data[j, i] > threshold_bonf) {
               threshold_bonf_count <- threshold_bonf_count + 1
               status <- 1
          }
     }
     conclude_skat_mx[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx[index, 2] <- threshold_bonf_count
     conclude_skat_mx[index, 3] <- threshold_bonf_norm_count
     conclude_skat_mx[index, 4] <- status
}

# The last row is show summation of fequency.
conclude_skat_mx[index+1, 1] <- "SUM"
conclude_skat_mx[index+1, 2] <- sum(as.numeric(conclude_skat_mx[, 2]))
conclude_skat_mx[index+1, 3] <- sum(as.numeric(conclude_skat_mx[, 3]))
conclude_skat_mx[index+1, 4] <- sum(as.numeric(conclude_skat_mx[, 4]))

# Set Header to result matrix
colnames(conclude_skat_mx) <- c("Rep",
                                "> Bonferroni",
                                "> Bonferroni normal",
                                "Status")

cat(sprintf("Writing result to %s\n", conclude_out_file))
write.csv(conclude_skat_mx, file = conclude_out_file, row.names = FALSE)


