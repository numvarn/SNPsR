# filename : count-additive-model-1dis-snp.R
# Edit by : Phisan Sookhee, 20 DEC 2016 11:50
#-------------------------------------------------------------------------------

setwd("~/ResearchCode/SNPsR")

conclude_single_file <- "result/single-snps-test-glm/conclusion-p-value/gene-effect-0.2onedis-AGaye/conclusion-100Rep.csv"
conclude_out_file_pos <- "result/single-snps-test-glm/conclusion-p-value/gene-effect-0.2onedis-AGaye/conclusion-100Rep-count-pos.csv"
conclude_out_file_neg <- "result/single-snps-test-glm/conclusion-p-value/gene-effect-0.2onedis-AGaye/conclusion-100Rep-count-neg.csv"
snpsName_path <- "input/gen-snps/NameONLY13479.csv"

# Read input data from CSV file15
cat(sprintf("Reding file %s\n", conclude_single_file))
conclude_glm_data <- read.csv(conclude_single_file, header = TRUE)

# Read input data from CSV files
cat(sprintf("Reding file %s\n", snpsName_path))
snpsName_data <- read.csv(snpsName_path, header = TRUE)

# Dedine disease SNPs for Aditive Model 
dis_snp_name1 <- "rs1861759"

snpsName <- snpsName_data[,1]
dis_snp_pos <- match(dis_snp_name1, snpsName)

col_number <- ncol(conclude_glm_data)
row_number <- nrow(conclude_glm_data)

# Cut only consider 2 SNP and store in matrix - spns_mx
conclude_glm_data_consider <- conclude_glm_data[dis_snp_pos, ]

threshold_bonf <- -log10(0.05 / row_number)
threshold_bonf_norm <- -log10(0.05)

# Create matrix for store result.
conclude_skat_mx <- matrix(0, col_number, 3)

# @Step 1 Count Positive
for (i in 2:ncol(conclude_glm_data_consider)) {
     index <- i - 1
     cat(sprintf("Processing Replicated #%s\n", index))
     status <- 0
     
     if(as.numeric(conclude_glm_data_consider[i]) > threshold_bonf) {
          status <- 1
     }
     

     conclude_skat_mx[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx[index, 2] <- as.numeric(conclude_glm_data_consider[i])
     conclude_skat_mx[index, 3] <- status
}

# The last row is show summation of fequency.
conclude_skat_mx[i, 1] <- "SUM"
conclude_skat_mx[i, 2] <- "-"
conclude_skat_mx[i, 3] <- sum(as.numeric(conclude_skat_mx[, 3]))

# Set Header to result matrix
colnames(conclude_skat_mx) <- c("Rep",
                                "Value",
                                "Status")

cat(sprintf("Writing result to %s\n", conclude_out_file_pos))
write.csv(conclude_skat_mx, file = conclude_out_file_pos, row.names = FALSE)

#-------------------------------------------------------------------------------
# @Step 2 Count Nagetive
conclude_glm_data_consider <- conclude_glm_data[-dis_snp_pos, ]

# Create matrix for store result.
conclude_skat_mx_negative <- matrix(0, ncol(conclude_glm_data_consider), 5)

for (i in 2:ncol(conclude_glm_data_consider)) {
     index <- i - 1
     status <- 0
     greater_bonf <- 0
     greater_bonf_norm <- 0
     
     greater_bonf_norm <- replicated <- conclude_glm_data_consider[, i]
     
     greater_bonf <- length(replicated[replicated>threshold_bonf])
     greater_bonf_norm <- length(replicated[replicated>threshold_bonf_norm])
     
     if(greater_bonf > 0) {
          status <- 1
     }
     
     conclude_skat_mx_negative[index, 1] <- paste("P", index, sep = "")
     conclude_skat_mx_negative[index, 2] <- greater_bonf
     conclude_skat_mx_negative[index, 3] <- greater_bonf_norm
     conclude_skat_mx_negative[index, 4] <- status
     conclude_skat_mx_negative[index, 5] <- greater_bonf / col_number
}

# The last row is show summation of fequency.
index <- index + 1
conclude_skat_mx_negative[index, 1] <- "SUM"
conclude_skat_mx_negative[index, 2] <- sum(as.numeric(conclude_skat_mx_negative[, 2]))
conclude_skat_mx_negative[index, 3] <- sum(as.numeric(conclude_skat_mx_negative[, 3]))
conclude_skat_mx_negative[index, 4] <- sum(as.numeric(conclude_skat_mx_negative[, 4]))
conclude_skat_mx_negative[index, 5] <- sum(as.numeric(conclude_skat_mx_negative[, 5]))

# Set Header to result matrix
colnames(conclude_skat_mx_negative) <- c("Rep",
                                "> Bonferroni",
                                "> Normal",
                                "Stauts",
                                "FT")

cat(sprintf("Writing result to %s\n", conclude_out_file_neg))
write.csv(conclude_skat_mx_negative, file = conclude_out_file_neg, row.names = FALSE)




