setwd("~/ResearchCode/SNPsR")

conclude_skat_file <- "result/single-snps-test-glm/conclusion-p-value/conclusion100Rep.csv"
conclude_tp_out_file <- "result/single-snps-test-glm/conclusion-p-value/count100Rep-TP.csv"
conclude_fp_out_file <- "result/single-snps-test-glm/conclusion-p-value/count1000Rep-FT.csv"
geneName_path <- "input/gen-snps/gene-name.csv"

# Read input data from CSV files
cat(sprintf("Reding file %s\n", conclude_skat_file))
conclude_skat_data <- read.csv(conclude_skat_file, header = TRUE)

# Read input data from CSV files
cat(sprintf("Reding file %s\n", geneName_path))
geneName_data <- read.csv(geneName_path, header = TRUE)

# Define disease SNPs for SKAT Model
dis_gene_name1 <- "RBFOX1"
dis_gene_name2 <- "NOD2"

geneName <- geneName_data[,1]
dis_gene_pos1 <- match(dis_gene_name1, geneName)
dis_gene_pos2 <- match(dis_gene_name2, geneName)

col_number <- ncol(conclude_skat_data)
row_number <- nrow(conclude_skat_data)

# Cut only consider 2 SNP and store in matrix - gene_mx
conclude_skat_data_consider <- conclude_skat_data[c(dis_gene_pos1, dis_gene_pos2), ]

threshold_bonf <- -log10(0.05/row_number)
# threshold_bonf_norm <- -log10(0.05)

# Create matrix for store result.
conclude_skat_mx <- matrix(0, col_number-7 ,6)

for (i in 9:ncol(conclude_skat_data_consider)) {
     index <- i - 8
     cat(sprintf("Processing Replicated #%s\n", index))
     threshold_bonf_count_1 <- 0
     threshold_bonf_count_2 <- 0
     status <- 0
     
     if(conclude_skat_data_consider[1, i] < threshold_bonf) {
          threshold_bonf_count_1 <- 1
     }
     
     if(conclude_skat_data_consider[2, i] > threshold_bonf) {
          threshold_bonf_count_2 <- 1
     }
     
     if (threshold_bonf_count_1 == 1 && threshold_bonf_count_2 == 1) {
          status <- 1
     }

     conclude_skat_mx[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx[index, 2] <- conclude_skat_data_consider[1, i]
     conclude_skat_mx[index, 3] <- threshold_bonf_count_1
     conclude_skat_mx[index, 4] <- conclude_skat_data_consider[2, i]
     conclude_skat_mx[index, 5] <- threshold_bonf_count_2
     conclude_skat_mx[index, 6] <- status
}

# The last row is show summation of frequency.
conclude_skat_mx[index+1, 1] <- "SUM"
conclude_skat_mx[index+1, 2] <- "-"
conclude_skat_mx[index+1, 3] <- sum(as.numeric(conclude_skat_mx[, 3]))
conclude_skat_mx[index+1, 4] <- "-"
conclude_skat_mx[index+1, 5] <- sum(as.numeric(conclude_skat_mx[, 5]))
conclude_skat_mx[index+1, 6] <- sum(as.numeric(conclude_skat_mx[, 6]))

# Set Header to result matrix
colnames(conclude_skat_mx)      <- c("Rep",
                                     "Gene1",
                                     "Dis1",
                                     "Gene2",
                                     "Dis2",
                                     "Status")

cat(sprintf("Writing result to %s\n", conclude_tp_out_file))
write.csv(conclude_skat_mx, file = conclude_tp_out_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# Path 2
# Count FP
conclude_fp <- conclude_skat_data[c(-dis_gene_pos1, -dis_gene_pos2), ]

# Create matrix for store result
conclude_fp_mx <- matrix(0, ncol(conclude_fp) - 7, 4)

for(i in 9:ncol(conclude_fp)) {
     index <- i - 8
     cat(sprintf("Processing Replicated #%s\n", index))
     
     fp_count <- 0
     status <- 0
     for(j in 1:nrow(conclude_fp)) {
          p_val <- as.numeric(conclude_fp[j, i])
          if (threshold_bonf < p_val) {
               fp_count <- fp_count + 1
          }
     }
     
     if(fp_count > 0) {
          status <- 1
     }
     
     conclude_fp_mx[index, ] <- c(paste("P", index), fp_count, status, fp_count/nrow(conclude_fp)) 
     
}

conclude_fp_mx[index+1, ] <- c("Sum",
                         sum(as.numeric(conclude_fp_mx[, 2])),
                         sum(as.numeric(conclude_fp_mx[, 3])),
                         sum(as.numeric(conclude_fp_mx[, 4]))
)

# Set Header to result matrix
colnames(conclude_fp_mx) <- c("Rep",
                              ">Bon",
                              "Status",
                              "FT")

cat(sprintf("Writing result to %s\n", conclude_fp_out_file))
write.csv(conclude_fp_mx, file = conclude_fp_out_file, row.names = FALSE)
