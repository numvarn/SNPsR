# Config values
setwd("~/ResearchCode/SNPsR")

# case_control_file <- commandArgs(TRUE)
case_control_file <- "result/gen-snps/case-control-replicated/0.0/1.csv"

#get replicated for CSV file
outfile_name <- basename(case_control_file)
outfile_p_value <- paste("result/single-snps-test/0.0/", outfile_name, sep = "")
     
replicated <- read.csv(case_control_file)

#get number of SNPs 
replicated_size <- ncol(replicated)
snps_number <- replicated_size - 3

#get Y
y.b <- replicated[, 1]

#create matrix for store result
p_mx <- matrix(0, ncol(replicated) - 3, 3)

#single SNPs test by General Linear Model (GLM)
snps_count <- 1
for (i in 4:replicated_size) {
     single_snps <- replicated[, i]
     
     #calculate p-value by using glm()
     results <- glm(y.b ~ single_snps, family = "binomial")
     deviance <- results$null.deviance - results$deviance
     df <- results$df.null - results$df.residual
     
     p_value <- 1 - pchisq(deviance, df)
     
     #store result into matrix
     p_mx[snps_count, 1] <- snps_count
     p_mx[snps_count, 2] <- p_value
     p_mx[snps_count, 3] <- -log10(p_value)
     
     snps_count <- snps_count + 1
}

# write results to CSV file
colnames(p_mx) <- c("SNPs NO.", "P-value", "-log10(P)")
write.csv(p_mx, file = outfile_p_value, row.names = FALSE)




