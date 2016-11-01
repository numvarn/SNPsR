# Config values
setwd("~/ResearchCode/SNPsR")
infile_result_clean <- "result/grouping-gene/05.GroupResultCleanManual.csv"
outfile_overlab_list <- "result/grouping-gene/06.OverlabSNPsList.csv"

overlab_list_mx <- matrix(NA, nrow = 1, ncol = 8)

cleaned_data <- read.csv(infile_result_clean, header =  TRUE)
for (i in 1:nrow(cleaned_data)) {
     if (!is.na(cleaned_data[i, 8])) {
          overlab_list_mx <- rbind(overlab_list_mx, sapply(cleaned_data[i, ], paste0, collapse=""))
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

# Write none-overlab file
colnames(overlab_list_mx) <- header
write.csv(overlab_list_mx[-1, ], file = outfile_overlab_list, row.names = FALSE)