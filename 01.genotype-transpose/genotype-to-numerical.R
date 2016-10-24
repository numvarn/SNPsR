setwd("~/ResearchCode/SNPsR")

geno_array_path <- "result/genotype-transpose/03.Geno_array_13479_Transpose_total.csv"
geno_data<- read.csv(geno_array_path, header = TRUE)

outfile <- "result/genotype-transpose/04.Geno_array_numerical_max.csv"

row_number <- nrow(geno_data)
col_number <- ncol(geno_data)

geno_numberical_mx <- matrix(0, row_number, col_number)

for (i in 1:col_number) {
     snp <- geno_data[, i]
     feq <- table(snp)
     a <- 0
     c <- 0
     g <- 0
     t <- 0
     
     if ("A" %in% snp) {
          a <- feq[["A"]][1]
     }
     
     if ("C" %in% snp) {
          c <- feq[["C"]][1]
     }
     
     if ("G" %in% snp) {
          g <- feq[["G"]][1]
     }
     
     if ("T" %in% snp) {
          t <- feq[["T"]][1]
     }
     
     #find index of max value
     ts <- c(a, c, g, t)
     max_index <- which(ts==max(ts))
     
     replace_str <- ""
     if (max_index == 1) {
          replace_str <- "A"
     } else if (max_index == 2) {
          replace_str <- "C"
     } else if (max_index == 3) {
          replace_str <- "G"
     } else if (max_index == 4) {
          replace_str <- "T"
     }
     
     cat(sprintf("Processing snps %s, %s\n", i, replace_str))
     
     #replace value in 0, 1
     for (j in 1:row_number) {
          if (geno_data[j, i] == replace_str) {
               geno_numberical_mx[j, i] <- 1
          } else {
               geno_numberical_mx[j, i] <- 0
          }
     }
}

# Write all data to CSV
write.csv(geno_numberical_mx, file = outfile, row.names = FALSE)

