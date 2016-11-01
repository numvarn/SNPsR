setwd("~/ResearchCode/SNPsR")

geno_array_path <- "input/genotype-transpose/Geno_array_13479_Transpose.csv"
geno_array_data<- read.csv(geno_array_path, header = FALSE)

snps_number <- ncol(geno_array_data)
rows_number <- nrow(geno_array_data)

for (i in 1:2) {
     snp <- geno_array_data[, i]
     feq = table(snp)
     
     double_caes <- c('AA', 'CC', 'GG', 'TT')
     
     v <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
     allele <- c('AA', 'AC', 'AG', 'AT', 
                 'CA', 'CC', 'CG', 'CT', 
                 'GA', 'GC', 'GG', 'GT', 
                 'TA', 'TC', 'TG', 'TT')
     names(v) <- allele
     
     for (key in allele) {
          if (key %in% snp) {
               v[key] <- feq[[key]][1]
          }
     }
     
     feq_summary <- v[v > 0]
     allele_summary <- c(names(feq_summary))
     
     n1 <- 0
     n2 <- 0
     n3 <- 0
     
     for (label in allele_summary) {
          if(label %in% double_caes && n1 == 0) {
               n1 = feq_summary[label]
          } else if(n2 == 0) {
               n2 = feq_summary[label]
          } else if(n2 != 0) {
               n3 = feq_summary[label]
          }
     }
   
     litter1 <- substr(allele_summary[2], 1, 1)
     litter2 <- substr(allele_summary[2], 2, 2)
     
     litter1_value = (2*n1) + n2
     litter2_value = (2*n3) + n2
     
     # replace litter with 1 by min value
     replace_litter <- ''
     if (litter1_value <= litter2_value) {
          replace_litter <- litter1
     } else {
          replace_litter <- litter2
     }
     
     # print value for dobule check
     print(feq_summary)
     cat(sprintf("n1=%s, n2=%s, n3=%s\n", n1, n2, n3))
     cat(sprintf("litter1: %s, litter2: %s\n", litter1, litter2))
     cat(sprintf("replace litter : %s\n", replace_litter))
     
}

