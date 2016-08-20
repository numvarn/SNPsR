# Clean Overlab Gene
# By unsing compairation between base-pair point and start, stop
# Edit by Phisan Sookkhe, 20 Aug 2016 - 12:15

# Config values
setwd("~/ResearchCode/SNPsR")

infile_result_overlab <- "result/grouping-gene/01.GroupResultOverlab.csv"
outfile_result_clean <- "result/grouping-gene/03.GroupResultClean.csv"
outfile_snp_on_gene <- "result/grouping-gene/04.SNPsonGeneClean.csv"

# Clean Overlab SNPs
result_mx <- read.csv(infile_result_overlab, header = TRUE)

result_mx_clean <- matrix("NA", nrow = 1, ncol = 8)
for (i in 1:nrow(result_mx)) {
     if (!is.na(result_mx[i, 8])) {
          crr_base_pair <- as.numeric(result_mx[i, 5]) 
          gene_1 <- result_mx[i-1, ]
          gene_2 <- result_mx[i, ]
          
          s1 <- as.numeric(gene_1[4])
          e1 <- as.numeric(gene_1[6])
          
          s2 <- as.numeric(gene_2[4])
          e2 <- as.numeric(gene_2[6])
          
          if (s1 < s2 && e2 < e1) {
               # @sub-set overlab case
               result_mx[i-1, 8] <- as.character("TRUE")
               result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i-1, ], paste0, collapse=""))
          } else {
               # @nornal overlab case
               
               # @Method - 1 : by using difference between start - stop
#                dt_bp_gene_1 <- abs(e1 - crr_base_pair)
#                dt_bp_gene_2 <- abs(s2 - crr_base_pair)
               
               # @Method - 2 : by using difference between median
               # @By experimental result : method - 2 is better than method - 1 
               dt_bp_gene_1 <- abs(median(s1:e1) - crr_base_pair)
               dt_bp_gene_2 <- abs(median(s2:e2) - crr_base_pair)
               
               if (dt_bp_gene_1 < dt_bp_gene_2) {
                    result_mx[i-1, 8] <- as.character("TRUE")
                    result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i-1, ], paste0, collapse=""))
               } 
               else {
                    result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i, ], paste0, collapse=""))
               }
          }
     } else {
          if (i < nrow(result_mx) && is.na(result_mx[i + 1, 8])) {
               result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i, ], paste0, collapse=""))
          }
     }
}

header <- c("no.",
            "rs number",
            "Bp old from Agaye", 
            "start", 
            "BP new", 
            "stop", 
            "symbols", 
            "overlab")

# Write none-overlab file
colnames(result_mx_clean) <- header
write.csv(result_mx_clean[-1, ], file = outfile_result_clean, row.names = FALSE)

# Count SNPs on earch Gene
cleaned <- read.csv(outfile_result_clean, header = TRUE)
snp_on_gene <- matrix(NA, 1, 3)
gene_count <- 0
for (i in 1:nrow(cleaned)) {
     symb <- cleaned[i, 7]
     if (!is.na(symb) && symb != "") {
          index <- 0
          index <- match(symb, as.vector(snp_on_gene[, 2]))
          if (is.na(index)) {
               gene_count <- gene_count + 1
               snp_on_gene <- rbind(snp_on_gene, c(gene_count, as.character(symb), 1))  
          } else {
               snp_on_gene[index, 3] <- as.numeric(snp_on_gene[index, 3]) + 1
          }
     }
}

# Write Result Overlab to CSV file ************************************
colnames(snp_on_gene) <- c("no.", "Gene Name", "SNPs Count")
write.csv(snp_on_gene[-1, ], file = outfile_snp_on_gene, row.names = FALSE)




