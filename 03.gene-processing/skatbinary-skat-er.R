# SKATBinary-SKAT-Hybrid (Defualt Parameter setting)
# *****************************************************************************#

library(SKAT)

# Config values
setwd("~/ResearchCode/SNPsR")
# case_control_file <- commandArgs(TRUE)
case_control_file <- "result/gen-snps/test/1.csv"

# check file-path from command line argument
#if (case_control_file[1] != "") {
if (case_control_file != "") {
     filename <- basename(case_control_file[1])
     start_stop_file <- "result/grouping-gene/10.GroupingComplete.csv"
     outfile_p_value <- paste("result/gen-snps/test/out/", filename, sep = "")
     
     z_temp <- read.csv(case_control_file, header = TRUE)
     start_stop = read.csv(start_stop_file, header = TRUE)
     yb <- z_temp[, 1]
     z <- z_temp[, 4:ncol(z_temp)]
     
     # delete error SNPs
     z <- z[, c(-5154)] #snps 13478
     
     # Convert Data Fram to matrix
     yb_mx <- as.matrix(yb)
     
     # Create matrix for store p-value result
     # pval <- matrix(NA, nrow(start_stop), 6)
     pval <- matrix(NA, 5, 6)
     
     # for(i in 1:nrow(start_stop)) {
     for(i in 1:5) {
          z_gene <- as.matrix(z[, start_stop[i, 3]:start_stop[i, 4]])
          obj <-
               SKAT_Null_Model(
                    yb_mx ~ 1,
                    out_type = "D"
               )
          p <- SKATBinary(
                    z_gene, 
                    obj, 
                    method.bin = "MA"
               )$p.value
          
          # cat(sprintf("Gene no.%s : %s , P-value = %s\n\n", i, start_stop[i, 2], p))
          pval[i, 1] <- i
          pval[i, 2] <- as.character(start_stop[i, 2])
          pval[i, 3] <- start_stop[i, 3]
          pval[i, 4] <- start_stop[i, 4]
          pval[i, 5] <- p
          pval[i, 6] <- -log10(p)
     }  
     colnames(pval) <- c("Gene No.", "Gene Name", "Start", "Stop", "P", "-log10(P)")
     write.csv(pval, file = outfile_p_value, row.names = FALSE)
} else {
     print("Please, Enter case-control file path !!")
}