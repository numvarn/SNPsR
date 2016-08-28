# Config values
setwd("~/ResearchCode/SNPsR")
# snp_pos_file <- "input/grouping-gene/SNP_Pos_13479_cut_2_snps.csv"

snp_pos_file <- "input/grouping-gene/NewSNPsBP_Sirikanlaya.csv"
start_stop_file <- "input/grouping-gene/StartStopPosition.csv"

outfile_result_overlab <- "result/grouping-gene/01.GroupResultOverlab.csv"
outfile_snp_on_gene_overlab <- "result/grouping-gene/02.SNPsonGeneOverlab.csv"

snp_pos_data <- read.csv(snp_pos_file, header = TRUE)
start_stop_data <- read.csv(start_stop_file, header = TRUE)

result <- c()
result_mx <- matrix(c("","","", "", "", "", "", ""), 1, 8)
snp_on_gene <- matrix(c("", "", ""), 1, 3)
gene_count <- 0

for (i in 1:nrow(snp_pos_data)) {
     cat(sprintf("Processing : %s - %s\n", i, snp_pos_data[i, 2]))
     found_count <- 0
     found <- FALSE
     
     base_pair <- snp_pos_data[i, 4]
     for (j in 1:nrow(start_stop_data)) {
          start_found <- 0
          stop_found <- 0
          symb <- NA
          overlab <- ""
          
          start_point <- start_stop_data[j, 1]
          stop_point <- start_stop_data[j, 2]
          
          if (start_point <= base_pair && base_pair <= stop_point) {
               found <- TRUE
               found_count <- found_count + 1
               
               start_found <- start_point
               stop_found <- stop_point
               symb <- start_stop_data[j, 3]
               
               # Overlab case
               if (found_count > 1) {
                    overlab <- "TRUE"
               }
               
               result <- c(snp_pos_data[i, 1],
                           as.character(snp_pos_data[i, 2]),
                           snp_pos_data[i, 3],
                           start_found,
                           snp_pos_data[i, 4], #base_pair
                           stop_found,
                           as.character(symb),
                           as.character(overlab))
               
               result_mx <- rbind(result_mx, result)
          }
          
          # Count SNP on each Gene
          if (!is.na(symb)) {
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
     
     if (found == FALSE) {
          result <- c(snp_pos_data[i, 1],
                      as.character(snp_pos_data[i, 2]),
                      snp_pos_data[i, 3],
                      as.character(""),
                      snp_pos_data[i, 4], #base_pair
                      as.character(""),
                      as.character(""),
                      as.character(""))
          
          result_mx <- rbind(result_mx, result)
     }
}

# Write Result Overlab to CSV file ************************************
colnames(snp_on_gene) <- c("no.", "Gene Name", "SNPs Count")
write.csv(snp_on_gene[-1, ], file = outfile_snp_on_gene_overlab, row.names = FALSE)

# Write overlab file
header <- c("no.",
            "rs number",
            "Bp old from Agaye", 
            "start", 
            "BP new", 
            "stop", 
            "symbols", 
            "overlab")

colnames(result_mx) <- header
write.csv(result_mx[-1, ], file = outfile_result_overlab, row.names = FALSE)




