# Code by Phisan Sookkhee and Sirikanlaya Sookkhee
# Edited 13 AUG 2016

# Config values
setwd("~/ResearchCode/SNPsR")
gene_pair_path <- "input/StartStopAnalysis.csv"
outfile_start_stop <- "result/StartStopResult.csv"

# Read file to data frame
gene_pair_data <- read.csv(gene_pair_path, header = TRUE)
number_row <- nrow(gene_pair_data)

current_gene <- ""
start_gene <- ""
start_point <- 0
stop_point <- 0
start_stop_vector <- c()

for (i in 1:number_row) {
     switch_flag = FALSE
     current_gene = as.character(gene_pair_data[i, 3])
     if (start_gene != current_gene) {
          # found stop point
          stop_point <- i
          switch_flag = TRUE
     } else if (start_gene == current_gene) {
          stop_point <- i
     } 
     
     if (switch_flag == TRUE) {
          start_stop_vector <- append(start_stop_vector, c(start_gene, start_point, stop_point))
          start_point <- i
          start_gene <- current_gene
     }
     
     if (i < 50) {
          cat(sprintf("%s : %s : %s, %s\n", i, as.character(gene_pair_data[i, 3]), start_point, stop_point)) 
     }
}

start_stop_matrix <- matrix(start_stop_vector, length(start_stop_vector)/3, 3, byrow = TRUE)

colnames(start_stop_matrix) <- c("Gene", "Start", "Stop")
write.csv(start_stop_matrix, file = outfile_start_stop, row.names = FALSE)


