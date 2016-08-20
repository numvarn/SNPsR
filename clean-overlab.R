# Config values
setwd("~/ResearchCode/SNPsR")

infile_result_overlab <- "result/GroupingGene/GroupResultOverlab.csv"
outfile_result_clean <- "result/GroupingGene/GroupResultClean.csv"
outfile_snp_on_gene <- "result/GroupingGene/SNPsonGeneClean.csv"

# Clean Overlab SNPs
result_mx <- read.csv(infile_result_overlab, header = TRUE)

result_mx_clean <- matrix("NA", nrow = 1, ncol = 8)
for (i in 1:nrow(result_mx)) {
     if (!is.na(result_mx[i, 8])) {
          crr_base_pair <- as.numeric(result_mx[i, 5]) 
          gene_1 <- result_mx[i-1, ]
          gene_2 <- result_mx[i, ]
          
          median_gene_1 <- median(c(as.numeric(gene_1[4]):as.numeric(gene_1[6])))
          median_gene_2 <- median(c(as.numeric(gene_2[4]):as.numeric(gene_2[6])))
          
          dt_bp_gene_1 <- abs(median_gene_1 - crr_base_pair)
          dt_bp_gene_2 <- abs(median_gene_2 - crr_base_pair)
          
          if (dt_bp_gene_1 < dt_bp_gene_2) {
               result_mx[i-1, 8] <- as.character("TRUE")
               result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i-1, ], paste0, collapse=""))
          } 
          else {
               result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i, ], paste0, collapse=""))
          }
     } else {
          if (i < nrow(result_mx) && is.na(result_mx[i + 1, 8])) {
               result_mx_clean <- rbind(result_mx_clean, sapply(result_mx[i, ], paste0, collapse=""))
          }
     }
}

# Write none-overlab file
header <- c("no.",
            "rs number",
            "Bp old from Agaye", 
            "start", 
            "BP new", 
            "stop", 
            "symbols", 
            "overlab")

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




