setwd("~/ResearchCode/SNPsR")

snps_name_file <- "input/check-snps/NameONLY13479.csv"
snps_basepair_file <- "input/check-snps/SNP_Pos_13479.csv"
outfile_result <- "result/check-snps/NewSNPsBP_Sirikanlaya.csv"

snps_name_data <- read.csv(snps_name_file, header = TRUE)
snps_basepair_data <- read.csv(snps_basepair_file, header = TRUE)

snps_main_data <- matrix(NA, nrow(snps_name_data), 3)
colnames(snps_main_data) <- c("No.", "SNPs", "BP")

for (i in 1:nrow(snps_name_data)) {
     bp <- NA
     current_snps <- as.character(snps_name_data[i, 1])
     for (j in 1:nrow(snps_basepair_data)) {
          check_snps <- as.character(snps_basepair_data[j, 2])
          if (current_snps == check_snps) {
               bp <- snps_basepair_data[j, 4]
               break
          }
     }
     snps_main_data[i, ] <- c(i, as.character(snps_name_data[i, 1]), bp)
}


write.csv(snps_main_data, file = outfile_result, row.names = FALSE)
