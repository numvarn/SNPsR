# For commbine two replicated.
# Original replicated is 1500 case and 1500 control.
# New replicated is 3000 case and 3000 control.
# Edited by Sirikanlaya Sookkhee, 25 May 2017 20:45 PM.

setwd("~/Desktop/concate-file")

outfile_path <- "conbined-replicated/" 


# Start program
# by get command line argument
# arg[1] is start number of replicated
# arg[2] is stop number of replicated
replicate_range <- commandArgs(TRUE)

if (length(replicate_range) == 2 && replicate_range[1] < replicate_range[2]) {
     cat(sprintf("start number : %s\n", replicate_range[1]))
     cat(sprintf("end number : %s\n", replicate_range[2]))


     start_filename <- as.numeric(replicate_range[1])
     end_filename <- as.numeric(replicate_range[2])
     
     #create odd number list
     odd_list <- seq(start_filename, end_filename, 2)
     
     for(i in odd_list) {
          filename_1 <- paste(i, ".csv", sep = "")
          filename_2 <- paste(i+1, ".csv", sep = "")
     
          original_rep_1 <- read.csv(filename_1, header = TRUE)
          original_rep_2 <- read.csv(filename_2, header = TRUE)
     
          new_rep <- rbind(original_rep_1, original_rep_2)
     
          outfile <- paste(outfile_path, i, ".csv", sep = "")
     
          # Write all data to CSV
          write.csv(new_rep, file = outfile, row.names = FALSE)
     }
}