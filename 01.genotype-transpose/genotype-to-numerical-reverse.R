setwd("~/ResearchCode/SNPsR")

geno_array_path <- "result/genotype-transpose/04.Geno_array_numerical_max.csv"
geno_data<- read.csv(geno_array_path, header = TRUE)

outfile <- "result/genotype-transpose/05.Geno_array_numerical_min.csv"

row_number <- nrow(geno_data)
col_number <- ncol(geno_data)

reverse_mx = matrix(0, row_number, col_number)

for (i in 1:col_number) {
     snp <- geno_data[, i]
     for (j in 1:row_number) {
          if (geno_data[j, i] == 1) {
               reverse_mx[j, i] <- 0
          } else {
               reverse_mx[j, i] <- 1
          }
     }
     
     cat(sprintf("Processing snps %s\n", i))
}

# Write all data to CSV
write.csv(reverse_mx, file = outfile, row.names = FALSE)

