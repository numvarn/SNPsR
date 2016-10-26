setwd("~/ResearchCode/SNPsR")

conclude_skat_file <- "result/skat/conclusion-p-value/conclusion-skat-1000.csv"
conclude_out_file <- "result/skat/conclusion-p-value/conclusion-skat1000_threshold.csv"

conclude_skat_data <- read.csv(conclude_skat_file, header = TRUE)
conclude_mx <- matrix(0, nrow(conclude_skat_data), 3)

col_number <- ncol(conclude_skat_data)
row_number <- nrow(conclude_skat_data)

threshold_bonf <- -log10(0.05 / row_number)
threshold_bonf_norm <- -log10(0.05)

for (i in 1:nrow(conclude_skat_data)) {
     row <- conclude_skat_data[i, 9:col_number]
     
     bonferroni_threshold <- row[row > threshold_bonf]
     bonferroni_norm_threshold <- row[row > threshold_bonf_norm]
     
     conclude_mx[i, 1] <- as.character(conclude_skat_data[i, 2])
     conclude_mx[i, 2] <- length(bonferroni_threshold)
     conclude_mx[i, 3] <- length(bonferroni_norm_threshold)
     
     cat(sprintf("%s : %s\n", i, conclude_mx[i, 1]))
}

colnames(conclude_mx) <- c("Gene Name", 
                           "Bonferroni threshold", 
                           "Bonferroni normal threshold")

write.csv(conclude_mx, file = conclude_out_file, row.names = FALSE)