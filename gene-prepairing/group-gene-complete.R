# Config values
setwd("~/ResearchCode/SNPsR")
group_gene_file <- "result/grouping-gene/05.GroupResultCleanManual.csv"
outfile_group_complete <- "result/grouping-gene/07.GroupingComplete.csv"

group_gene_data <- read.csv(group_gene_file, header = TRUE )
group_gene_complete_mx <- matrix(NA, nrow = nrow(group_gene_data), ncol = 8)
symb <- c()

# @Step 1 : self-declare gene name for NA value
gene_na_count <- 0
i <- 1
while (i <= nrow(group_gene_data)) {
     if (as.character(group_gene_data[i, 7]) == "") {
          gene_na_count <- gene_na_count + 1
          gene_name <- paste("IG", gene_na_count, sep = "")

          for (x in 1:ncol(group_gene_data)) {
               group_gene_complete_mx[i, x] <- as.character(group_gene_data[i, x])
          }
          group_gene_complete_mx[i, 7] <- gene_name
          
          j <- i + 1
          while (as.character(group_gene_data[j, 7]) == "") {
               for (x in 1:ncol(group_gene_data)) {
                    group_gene_complete_mx[j, x] <- as.character(group_gene_data[j, x])
               }
               group_gene_complete_mx[j, 7] <- gene_name
               j <- j + 1
          }
          i <- j
     }else {
          for (x in 1:ncol(group_gene_data)) {
               group_gene_complete_mx[i, x] <- as.character(group_gene_data[i, x])
          }
          i <- i + 1
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

# Write Result to CSV file
colnames(group_gene_complete_mx) <- header
write.csv(group_gene_complete_mx, file = outfile_group_complete, row.names = FALSE)
