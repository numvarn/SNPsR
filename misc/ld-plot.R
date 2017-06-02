# @Section 1 : Configure values
setwd("~/ResearchCode/SNPsR")

ld_value_path <- "/Volumes/Sirikanlaya/DataFromAGaye/LD-Split/500/2-rs1861759.csv"

cat(sprintf("Reading file LD Value : %s", ld_value_path))
ld_value_data <- read.csv(ld_value_path, header = TRUE)

result_mx <- matrix(0, nrow(ld_value_data), 2)
for (i in 1:nrow(ld_value_data)) {
     result_mx[i, 1] <- i
     result_mx[i, 2] <- ld_value_data[i, 1]
}

# Plot Graph
library(ggplot2)
dat <- data.frame(
          no = result_mx[, 1],
          ld = result_mx[, 2]
     )

# Plot without Lable
ggplot(dat, aes(x=no, y=ld)) + 
     geom_point(shape=19, alpha=1) +
     geom_hline(yintercept=0.5, linetype="dashed", color = "red") + 
     geom_line() #+
     #geom_text(aes(label=no), size=3, hjust=0, vjust=2, check_overlap=TRUE)

# Plot with Lable
# ggplot(dat, aes(x=no, y=ld)) + geom_point(shape=16) + 
#      geom_text(aes(label=no), size=1, hjust=0, vjust=2, check_overlap=TRUE)





     
     