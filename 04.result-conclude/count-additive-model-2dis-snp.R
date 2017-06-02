setwd("~/ResearchCode/SNPsR")

conclude_file <- "result/single-snps-test-glm/conclusion-p-value/conclusion1000Rep.csv"
conclude_out_file <- "result/single-snps-test-glm/conclusion-p-value/conclusion1000Rep-count.csv"
conclude_fp_out_file <- "result/single-snps-test-glm/conclusion-p-value/conclusion1000Rep-fp.csv"
snpsName_path <- "input/gen-snps/NameONLY13479.csv"

# Read input data from CSV files
cat(sprintf("Reding file %s\n", conclude_file))
conclude_glm_data <- read.csv(conclude_file, header = TRUE)

# Read input data from CSV files
cat(sprintf("Reding file %s\n", snpsName_path))
snpsName_data <- read.csv(snpsName_path, header = TRUE)

# Dedine disease SNPs for Aditive Model
dis_snp_name1 <- "rs3789038"
dis_snp_name2 <- "rs3785142"

snpsName <- snpsName_data[,1]
dis_snp_pos1 <- match(dis_snp_name1, snpsName)
dis_snp_pos2 <- match(dis_snp_name2, snpsName)

col_number <- ncol(conclude_glm_data)
row_number <- nrow(conclude_glm_data)

# Cut only consider 2 SNP and store in matrix - spns_mx
conclude_glm_data_consider <- conclude_glm_data[c(dis_snp_pos1, dis_snp_pos2), ]

threshold_bonf <- -log10(0.05 / row_number)
# threshold_bonf_norm <- -log10(0.05)

# Create matrix for store result.
conclude_skat_mx <- matrix(0, col_number, 6)

for (i in 2:ncol(conclude_glm_data_consider)) {
     index <- i - 1
     cat(sprintf("Processing Replicated #%s\n", index))
     threshold_bonf_count_1 <- 0
     threshold_bonf_count_2 <- 0
     status <- 0
     
     if(conclude_glm_data_consider[1, i] > threshold_bonf) {
          threshold_bonf_count_1 <- 1
     }
     
     if(conclude_glm_data_consider[2, i] > threshold_bonf) {
          threshold_bonf_count_2 <- 1
     }
     
     if (threshold_bonf_count_1 == 0 && threshold_bonf_count_2 == 0) {
          status <- 0
     } else {
          status <- 1
     }

     conclude_skat_mx[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx[index, 2] <- conclude_glm_data_consider[1, i]
     conclude_skat_mx[index, 3] <- threshold_bonf_count_1
     conclude_skat_mx[index, 4] <- conclude_glm_data_consider[2, i]
     conclude_skat_mx[index, 5] <- threshold_bonf_count_2
     conclude_skat_mx[index, 6] <- status
}

# The last row is show summation of fequency.
conclude_skat_mx[index+1, 1] <- "SUM"
conclude_skat_mx[index+1, 2] <- "-"
conclude_skat_mx[index+1, 3] <- sum(as.numeric(conclude_skat_mx[, 3]))
conclude_skat_mx[index+1, 4] <- "-"
conclude_skat_mx[index+1, 5] <- sum(as.numeric(conclude_skat_mx[, 5]))
conclude_skat_mx[index+1, 6] <- sum(as.numeric(conclude_skat_mx[, 6]))

# Set Header to result matrix
colnames(conclude_skat_mx) <- c("Rep",
                                "Value1",
                                "Dis1",
                                "Value2",
                                "Dis2",
                                "Status")

cat(sprintf("Writing result to %s\n", conclude_out_file))
write.csv(conclude_skat_mx, file = conclude_out_file, row.names = FALSE)


# ------------------------------------------------------------------------------
# Path 2
# Count FP
conclude_glm_data_fp <- conclude_glm_data[c(-dis_snp_pos1, -dis_snp_pos2), ]

conclude_tp_mx <- matrix(0, ncol(conclude_glm_data_fp), 4)

for(i in 2:ncol(conclude_glm_data_fp)) {
     fp_count <- 0
     status <- 0
     for(j in 1:nrow(conclude_glm_data_fp)) {
          p_val <- as.numeric(conclude_glm_data_fp[j, i])
          if (threshold_bonf < p_val) {
               fp_count <- fp_count + 1
          }
     }
     
     if(fp_count > 0) {
          status <- 1
     }
     
     conclude_tp_mx[i-1, ] <- c(paste("P", i-1), fp_count, status, fp_count/nrow(conclude_glm_data_fp)) 
     
}

conclude_tp_mx[i, ] <- c("Sum",
                           sum(as.numeric(conclude_tp_mx[, 2])),
                           sum(as.numeric(conclude_tp_mx[, 3])),
                           sum(as.numeric(conclude_tp_mx[, 4]))
                         )

# Set Header to result matrix
colnames(conclude_tp_mx) <- c("Rep",
                                ">Bon",
                                "Status",
                                "FT")

cat(sprintf("Writing result to %s\n", conclude_fp_out_file))
write.csv(conclude_tp_mx, file = conclude_fp_out_file, row.names = FALSE)








