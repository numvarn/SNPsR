library(SKAT)

# Config values
setwd("~/ResearchCode/SNPsR")
z_file <- "input/skat/Z.csv"
yb_file <- "input/skat/yb.csv"
start_stop_file <- "input/skat/StartStopResult.csv"
outfile_p_value <- "result/skat/pvalue.csv"

z <- read.csv(z_file, header = FALSE)
yb <- read.csv(yb_file, header = FALSE)
start_stop = read.csv(start_stop_file, header = TRUE)

# Convert Data Fram to matrix
yb_mx <- as.matrix(yb[, 1])
pval <- matrix(NA, nrow(start_stop), 6)

for(i in 1:nrow(start_stop)) {
#for(i in 1:5) {
     z_gene = as.matrix(z[, start_stop[i, 3]:start_stop[i, 4]])
     obj <-
          SKAT_Null_Model(
               yb_mx ~ 1,
               out_type = "D",
               n.Resampling = 2000,
               type.Resampling = "bootstrap"
          )
     p <- SKAT(z_gene, obj, kernel = "linear.weighted")$p.value
     
     cat(sprintf("Gene no.%s : %s , P-value = %s\n\n", i, start_stop[i, 2], p))
     
     pval[i, 1] <- i
     pval[i, 2] <- as.character(start_stop[i, 2])
     pval[i, 3] <- start_stop[i, 3]
     pval[i, 4] <- start_stop[i, 4]
     pval[i, 5] <- p
     pval[i, 6] <- -log10(p)
}

colnames(pval) <- c("Gene No.", "Gene Name", "Start", "Stop", "P", "-log10(P)")
write.csv(pval, file = outfile_p_value, row.names = FALSE)

