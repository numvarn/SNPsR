#-------------------------------------------------------------------------------
#
#
#
#
#
#-------------------------------------------------------------------------------

# Configure value
setwd("~/ResearchCode/SNPsR")

geno_array_path <- "input/genotype-transpose/Geno_array_13479_Transpose.csv"

cat(sprintf("Read file %s\n", geno_array_path))
geno_array_data<- read.csv(geno_array_path, header = FALSE)

outfile1 <- "result/genotype-transpose/01.Geno_array_13479_Transpose_binary.csv"
outfile2 <- "result/genotype-transpose/02.Geno_array_13479_Transpose_char.csv"

snps_number <- ncol(geno_array_data)
rows_number <- nrow(geno_array_data)

transpose_mx <- matrix(0, rows_number * 2, snps_number)
transpose_mx_char <- matrix(0, rows_number * 2, snps_number)

# for (i in 123:123) {
for (i in 1:snps_number) {
     snp <- geno_array_data[, i]
     
     # Get fequency of allele in each SNPs
     feq = table(snp)
     
     double_caes <- c('AA', 'CC', 'GG', 'TT')
     
     v <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
     allele <- c('AA', 'AC', 'AG', 'AT', 
                 'CA', 'CC', 'CG', 'CT', 
                 'GA', 'GC', 'GG', 'GT', 
                 'TA', 'TC', 'TG', 'TT')
     names(v) <- allele
     
     # Store fequency value in vector
     for (key in allele) {
          if (key %in% snp) {
               v[key] <- feq[[key]][1]
          }
     }
     
     # Filter allele in SNPs by using value > 0
     feq_summary <- v[v > 0]
     allele_summary <- c(names(feq_summary))
     
     n1 <- 0
     n2 <- 0
     n3 <- 0
     
     # Find min, max value of allele
     # For case example. AA AG GG
     if (length(feq_summary) == 3) {
          alphabet1 <- substr(allele_summary[2], 1, 1)
          alphabet2 <- substr(allele_summary[2], 2, 2)
          
          n1 = as.integer(feq_summary[1])
          n2 = as.integer(feq_summary[2])
          n3 = as.integer(feq_summary[3])
          
     } else if (length(feq_summary) == 2) {
          # For case example. AA AG
          if (allele_summary[1] %in% double_caes && 
              !(allele_summary[2] %in% double_caes)) {
               alphabet1 <- substr(allele_summary[2], 1, 1)
               alphabet2 <- substr(allele_summary[2], 2, 2)
               
               n1 = as.integer(feq_summary[1])
          
          } # For case example. AG GG 
          else if (!(allele_summary[1] %in% double_caes) && 
                     allele_summary[2] %in% double_caes) {
               alphabet1 <- substr(allele_summary[1], 1, 1)
               alphabet2 <- substr(allele_summary[1], 2, 2)
               
               n2 = as.integer(feq_summary[1])

          } # For case example. AA GG 
          else if (allele_summary[1] %in% double_caes && 
                   allele_summary[2] %in% double_caes) {
               alphabet1 <- substr(allele_summary[1], 1, 1)
               alphabet2 <- substr(allele_summary[2], 1, 1)
               
               n1 = as.integer(feq_summary[1])

          }
          
          if (allele_summary[2] %in% double_caes) {
               n3 = as.integer(feq_summary[2])    
          } else {
               n2 = as.integer(feq_summary[2])
          }
     }
     
     # Formular for find min, max
     alphabet1_value = (2*n1) + n2
     alphabet2_value = (2*n3) + n2
     
     # replace litter with 1 by min value
     replace_alphabet <- ''
     if (alphabet1_value <= alphabet2_value) {
          replace_alphabet <- alphabet1
     } else {
          replace_alphabet <- alphabet2
     }
     
     # Final process
     # Split allele and replace to binary
     index <- 1
     for (allele_str in snp) {
          left_binary <- 0
          right_binary <- 0
          
          left <- substr(allele_str, 1, 1)
          right <- substr(allele_str, 2, 2)
          
          if (left == replace_alphabet) {
              left_binary <- 1
          }
          if (right == replace_alphabet) {
              right_binary <- 1
          }
          
          # Store result to transpose matrix
          # Process left
          transpose_mx[index, i] <- left_binary
          transpose_mx_char[index, i] <- left
          index <- index + 1
          
          # Process right
          transpose_mx[index, i] <- right_binary
          transpose_mx_char[index, i] <- right
          index <- index + 1
     }
     
     cat(sprintf("Processing SNPs #%s\n", i))
     
#      cat(sprintf("%s, %s, %s\n", n1, n2, n3))
#      print(feq_summary)
#      cat(sprintf("alphabet %s, %s\n", alphabet1, alphabet2))
#      cat(sprintf("replace : %s\n", replace_alphabet))
}

# Write all data to CSV
print("Write result to 01.Geno_array_13479_Transpose_binary.csv")
write.csv(transpose_mx, file = outfile1, row.names = FALSE)

print("Write result to 02.Geno_array_13479_Transpose_char.csv")
write.csv(transpose_mx_char, file = outfile2, row.names = FALSE)

