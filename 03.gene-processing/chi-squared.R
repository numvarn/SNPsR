# Config values
setwd("~/ResearchCode/SNPsR")

replicated_file <- commandArgs(TRUE)
filename <- basename(replicated_file[1])

outfile_p_value <- paste("result/chi-squared/0.0/", filename, sep = "")

replicated_data <- read.csv(replicated_file, header = TRUE)

y <- replicated_data[, 1]
p_mx <- matrix(0, ncol(replicated_data) - 3, 3)
count_snps <- 1

for (i in 4:ncol(replicated_data)) {
     snps <- replicated_data[, i]
     tbl <- table(y, snps)

     p_mx[count_snps, 1] <- count_snps
     p_mx[count_snps, 2] <- chisq.test(tbl)$p.value 
     p_mx[count_snps, 3] <- -log10(p_mx[count_snps, 2])
     
     count_snps <- count_snps + 1
}

colnames(p_mx) <- c("SNPs NO.", "P-value", "-log10(P)")
write.csv(p_mx, file = outfile_p_value, row.names = FALSE)


