library(SKAT)

# Config values
setwd("~/ResearchCode/SNPsR")

case_control_file <- "result/gen-snps/test/1.csv"
start_stop_file <- "result/grouping-gene/10.GroupingComplete.csv"

z_temp <- read.csv(case_control_file, header = TRUE)
start_stop = read.csv(start_stop_file, header = TRUE)

yb <- z_temp[, 1]
z <- z_temp[, 4:ncol(z_temp)]

# delete error SNPs
z <- z[, c(-5154)] #snps 13478

i <- 2
z_gene <- as.matrix(z[, start_stop[i, 3]:start_stop[i, 4]])

obj<-SKAT_Null_Model(yb ~ 1, out_type="D")

SKATBinary(z_gene, obj, method.bin="QA")$p.value