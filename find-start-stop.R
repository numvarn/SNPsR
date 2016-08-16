# Code by Phisan Sookkhee and Sirikanlaya Sookkhee
# Edited 13 AUG 2016 - 1.30 pm

# Config values
setwd("~/ResearchCode/SNPsR")
gene_pair_path <- "input/StartStopAnalysis.csv"
outfile_start_stop <- "result/StartStopResult.csv"

# Read file to data frame
gene_pair_data <- read.csv(gene_pair_path, header = TRUE)
number_row <- nrow(gene_pair_data)

start_point <- 1
stop_point <- 1
start_stop_vector <- c()
gene_number <- 0

for (i in 1:number_row) {
     current_gene <- as.character(gene_pair_data[i, 3])
     next_gene <- as.character(gene_pair_data[(i+1), 3])
     
     if (current_gene != next_gene && i < number_row) {
          gene_number <- gene_number + 1
          stop_point <- i
          start_stop_vector <- append(start_stop_vector, c(gene_number, as.character(gene_pair_data[i, 3]), start_point, stop_point))
          start_point <- i + 1
     }
}

start_stop_matrix <- matrix(start_stop_vector, length(start_stop_vector)/4, 4, byrow = TRUE)

colnames(start_stop_matrix) <- c("Gene No.", "Gene Name", "Start", "Stop")
write.csv(start_stop_matrix, file = outfile_start_stop, row.names = FALSE)


