# Auther Phisan Sookkhee

setwd("~/ResearchCode/SNPsR")

conclude_skat_file <- "result/skat/conclusion-p-value/conclusion-0.2-100Rep.csv"
conclude_out_file <- "result/skat/conclusion-p-value/conclusion-0.2-100Rep-out.csv"
geneName_path <- "input/gen-snps/gene-name.csv"

# Read input data from CSV files
cat(sprintf("Reding file %s\n", conclude_skat_file))
conclude_glm_data <- read.csv(conclude_skat_file, header = TRUE)

# Read input data from CSV files
cat(sprintf("Reding file %s\n", geneName_path))
snpsName_data <- read.csv(geneName_path, header = TRUE)

# Dedine disease SNPs for Aditive Model
dis_snp_name1 <- "RBFOX1"
dis_snp_name2 <- "NOD2"

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
conclude_skat_mx <- matrix(0, col_number-7, 6)

for (i in 9:ncol(conclude_glm_data_consider)) {
     index <- i - 8
     cat(sprintf("Processing Replicated #%s\n", index))
     threshold_bonf_count_1 <- 0
     threshold_bonf_count_2 <- 0
     status <- 0
     
     if(conclude_glm_data_consider[1, i] < threshold_bonf) {
          threshold_bonf_count_1 <- 1
     }
     
     if(conclude_glm_data_consider[2, i] > threshold_bonf) {
          threshold_bonf_count_2 <- 1
     }
     
     if (threshold_bonf_count_1 == 1 && threshold_bonf_count_2 == 1) {
          status <- 1
     }

     conclude_skat_mx[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx[index, 2] <- conclude_glm_data_consider[1, i]
     conclude_skat_mx[index, 3] <- conclude_glm_data_consider[2, i]
     conclude_skat_mx[index, 4] <- threshold_bonf_count_1
     conclude_skat_mx[index, 5] <- threshold_bonf_count_2
     conclude_skat_mx[index, 6] <- status
}

# The last row is show summation of fequency.
conclude_skat_mx[index+1, 1] <- "SUM"
conclude_skat_mx[index+1, 2] <- sum(as.numeric(conclude_skat_mx[, 2]))
conclude_skat_mx[index+1, 3] <- sum(as.numeric(conclude_skat_mx[, 3]))
conclude_skat_mx[index+1, 4] <- sum(as.numeric(conclude_skat_mx[, 4]))

# Set Header to result matrix
colnames(conclude_skat_mx) <- c("Rep",
                                "Gene1",
                                "Gene2",
                                "Dis1",
                                "Dis2",
                                "Status")

cat(sprintf("Writing result to %s\n", conclude_out_file))
write.csv(conclude_skat_mx, file = conclude_out_file, row.names = FALSE)


