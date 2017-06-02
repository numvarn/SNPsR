# @Section 1 : Configure values
setwd("~/ResearchCode/SNPsR")

dis_snp <- 'rs1861759'

ld_value_path <- "/Volumes/Sirikanlaya/DataFromAGaye/LDData-new.csv"
snpsName_path <- "input/gen-snps/NameONLY13479.csv"

# Read data from CSV
cat(sprintf("Reading file LD Value : %s", ld_value_path))
ld_value_data <- read.csv(ld_value_path, header = TRUE)

cat(sprintf("Reading file SNPs Name file : %s", snpsName_path))
snpsName_data <- read.csv(snpsName_path, header = TRUE)
snpsName_data <- snpsName_data[, 1]

dis_snp_pos <- match(dis_snp, snpsName_data)

selected_snps <- ld_value_data[, dis_snp_pos]

max_length <- length(selected_snps)
start_index <- 1
stop_index <- 1
filename <- 1

for (index in 1:max_length) {
     if (index %% 500 == 0 || index == max_length) {
          stop_index <- index
          cat(sprintf("file %s.csv\tsatrt : %s , stop : %s\n", filename, start_index, stop_index))
          
          result_mx <- as.matrix(selected_snps[start_index:stop_index])
          colnames(result_mx) <- c(dis_snp)
          
          outfile <- paste("/Volumes/Sirikanlaya/DataFromAGaye/LD-Split/", filename, "-", dis_snp, ".csv", sep = "")
          
          # Write all data to CSV
          write.csv(result_mx, file = outfile, row.names = FALSE)
          
          start_index <- stop_index + 1
          filename <- filename + 1
     }
}
