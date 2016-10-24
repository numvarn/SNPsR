setwd("~/ResearchCode/SNPsR")

geno_array_path <- "input/genotype-transpose/Geno_array_13479_Transpose.csv"
geno_array_data<- read.csv(geno_array_path, header = FALSE)

outfile_left <- "result/genotype-transpose/01.Geno_array_13479_Transpose_left.csv"
outfile_right <- "result/genotype-transpose/02.Geno_array_13479_Transpose_right.csv"
outfile_total <- "result/genotype-transpose/03.Geno_array_13479_Transpose_total.csv"

snps <- ncol(geno_array_data)
rows <- nrow(geno_array_data)

geno_left_mx <- matrix(NA, rows, snps)
geno_right_mx <- matrix(NA, rows, snps)

for (i in 1:rows) {
     line <- c()
     for (j in 1:snps) {
          char <- as.character(geno_array_data[i, j])
          line <- c(line, char)
     }
     geno_left_mx[i, ] <- substr(line, 1, 1)
     geno_right_mx[i, ] <- substr(line, 2, 2)
     
     cat(sprintf("Processing row : %s\n", i))
}

total_mx = rbind(geno_left_mx, geno_right_mx)

# Write all data to CSV
write.csv(geno_left_mx, file = outfile_left, row.names = FALSE)
write.csv(geno_right_mx, file = outfile_right, row.names = FALSE)
write.csv(total_mx, file = outfile_total, row.names = FALSE)


