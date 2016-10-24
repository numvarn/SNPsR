# Config values
setwd("~/ResearchCode/SNPsR")
group_gene_file <- "result/grouping-gene/05.GroupResultCleanManual.csv"

outfile_group_complete <- "result/grouping-gene/07.GroupingComplete.csv"
outfile_group_start_stop <- "result/grouping-gene/08.GroupingStartStop.csv"
outfile_grouping_new <- "result/grouping-gene/09.GroupingNewSingle.csv"
outfile_grouping_conplete <- "result/grouping-gene/10.GroupingComplete.csv"

group_gene_data <- read.csv(group_gene_file, header = TRUE )
group_gene_mx <- as.matrix(group_gene_data)

# ******************************************************************************
# @Step 1 : self-declare gene name for NA value
gene_na_count <- 0
i <- 1

while (i <= nrow(group_gene_data)) {
     if (as.character(group_gene_data[i, 7]) == "") {
          gene_na_count <- gene_na_count + 1
          gene_name <- paste("IG", gene_na_count, sep = "")
          group_gene_mx[i, 7] <- gene_name
          
          j <- i + 1
          while (as.character(group_gene_data[j, 7]) == "") {
               group_gene_mx[j, 7] <- gene_name
               j <- j + 1
          }
          i <- j
     }else {
          i <- i + 1
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

# Write Result to CSV file
colnames(group_gene_mx) <- header
write.csv(group_gene_mx, file = outfile_group_complete, row.names = FALSE)

# ******************************************************************************
# @Step 2 : Grouping Gene
#           by find start-stop of gene
number_row <- nrow(group_gene_mx)
start_point <- 1
stop_point <- 1
start_stop_vector <- c()
gene_number <- 0
start_stop_matrix <- matrix(NA, 1, 6)

for (i in 1:number_row) {
     current_gene <- as.character(group_gene_mx[i, 7])
     
     if (i+1 <= number_row) {
          next_gene <- as.character(group_gene_mx[(i+1), 7])
     } else {
          next_gene <- "LAST-GENE"
     }
     
     if (current_gene != next_gene) {
          gene_number <- gene_number + 1
          stop_point <- i
          
          start_bp <- group_gene_mx[start_point, 5]
          stop_bp <- group_gene_mx[stop_point, 5]
          
          # Check for IG case
          if (is.na(start_bp) && is.na(stop_bp)) {
               start_bp <- group_gene_mx[start_point, 5]
               stop_bp <- group_gene_mx[stop_point, 5]
          }
          
          # store result in matrix
          start_stop_matrix <- rbind(start_stop_matrix, c(gene_number, 
                                                          as.character(group_gene_mx[i, 7]), 
                                                          start_point, 
                                                          stop_point,
                                                          start_bp,
                                                          stop_bp))
          start_point <- i + 1
     } 
}

# Write Result to CSV file
colnames(start_stop_matrix) <- c("Group No.", 
                                 "Group Name", 
                                 "Start No.", 
                                 "Stop No.",
                                 "Start BP",
                                 "Stop BP")

write.csv(start_stop_matrix[-1, ], file = outfile_group_start_stop, row.names = FALSE)


# ******************************************************************************
# @Step 3 : Grouping Continue Multiple Single SNPs on Single Gene to one group
number_row <- nrow(start_stop_matrix[-1, ])
grouping_mx <- start_stop_matrix[-1, ]
new_grouping_mx <- matrix(NA, 1, 7)
single_gene_group_count <- 1

i <- 1
while (i <= number_row) {
     # Find Single SNPs on Single Gene
     start_point <- grouping_mx[i, 3]
     stop_point <- grouping_mx[i, 4]
     members <- c()
     
     if (start_point != stop_point) {
          new_grouping_mx <- rbind(new_grouping_mx, grouping_mx[i, ])
          new_group_no <- nrow(new_grouping_mx)
          new_grouping_mx[new_group_no, 1] <- new_group_no - 1
          new_grouping_mx[new_group_no, 7] <- ""
     } else {
          check <- TRUE
          j <- i
          single_gene_count <- 1
          start_bp <- grouping_mx[i, 5]
          stop_bp <- grouping_mx[i, 6]
          
          while (check) {
               curr_gene_name <- as.character(grouping_mx[j, 2])
               next_gene_name <- as.character(grouping_mx[j+1, 2])
               
               members <- c(members, curr_gene_name)
               
               next_gene_start <- grouping_mx[j+1, 3]
               next_gene_stop <- grouping_mx[j+1, 4]
               
               if (next_gene_start == next_gene_stop) {
                    stop_point <- next_gene_stop
                    j <- j + 1
                    single_gene_count <- single_gene_count + 1
                    stop_bp <- grouping_mx[j, 6]
                    if (!is.integer(match(next_gene_name, members))) {
                         members <- c(members, next_gene_name)
                    }
               } else {
                    i <- j
                    check <- FALSE
               }
          }
          
          if (single_gene_count > 1) {
               group_name <- paste("SingleGeneGroup", single_gene_group_count, sep = "")
               single_gene_group_count <- single_gene_group_count + 1
          } else {
               group_name <- grouping_mx[i, 2]
          }
          
          new_grouping_mx <- rbind(new_grouping_mx, c("", 
                                                      group_name, 
                                                      start_point, 
                                                      stop_point, 
                                                      start_bp, 
                                                      stop_bp,
                                                      toString(members)))
          new_group_no <- nrow(new_grouping_mx)
          new_grouping_mx[new_group_no, 1] <- new_group_no - 1
     }
     
     i <- i + 1
}

# Write Result to CSV file
colnames(new_grouping_mx) <- c("Group No.", 
                                 "Group Name", 
                                 "Start No.", 
                                 "Stop No.",
                                 "Start BP",
                                 "Stop BP", 
                                 "New Members")
write.csv(new_grouping_mx[-1, ], file = outfile_grouping_new, row.names = FALSE)

# ******************************************************************************
# @Step 4 : Grouping Gene by considering only single gene between two groups

grouping_mx_s4 <- new_grouping_mx[-1, ]
grouping_mx_complete <- matrix(NA, 1, 8)
number_row <- nrow(grouping_mx_s4)
i <- 1

for (i in 1:number_row) {
     # Find Single SNPs on Single Gene
     current_gene <- grouping_mx_s4[i, 2]
     start_point <- grouping_mx_s4[i, 3]
     stop_point <- grouping_mx_s4[i, 4]
     
     if (start_point == stop_point) {
          current_bp <- grouping_mx_s4[i, 5]
          left_stop_bp <- grouping_mx_s4[i-1, 6]
          right_strt_bp <- grouping_mx_s4[i+1, 5]
          
          dt_from_left <- as.numeric(current_bp) - as.numeric(left_stop_bp)
          dt_from_right <- as.numeric(right_strt_bp) - as.numeric(current_bp)
          
          if (dt_from_left <= dt_from_right) {
               index <- i - 1
               grouping_mx_s4[i-1, 6] <- current_bp
               grouping_mx_s4[i-1, 4] <- stop_point
          } else {
               index <- i + 1
               grouping_mx_s4[i+1, 5] <- current_bp
               grouping_mx_s4[i+1, 3] <- start_point
          }
          
          if (grouping_mx_s4[index, 7] == "") {
               grouping_mx_s4[index, 7] <- NA
          }
          
          member <- c(grouping_mx_s4[index, 2], 
                      grouping_mx_s4[index, 7], 
                      grouping_mx_s4[i, 7])
          member <- na.omit(member)
          grouping_mx_s4[index, 7] <- toString(member)
     }
}

# Clean single gene and run new number
for (i in 1:number_row) {
     start_point <- grouping_mx_s4[i, 3]
     stop_point <- grouping_mx_s4[i, 4]
     if (start_point != stop_point) {
          grouping_mx_complete <- rbind(grouping_mx_complete, grouping_mx_s4[i, ])
     }
}

for (i in 1:nrow(grouping_mx_complete)) {
     grouping_mx_complete[i, 1] <- (i - 1)
}

# Find median of start - stop point
number_row <- nrow(grouping_mx_complete)
for (i in 2:number_row) {
     start_bp <- grouping_mx_complete[i, 5]
     stop_bp <- grouping_mx_complete[i, 6]
     md <- median(start_bp:stop_bp)
     grouping_mx_complete[i, 8] <- md
}


# Write Result to CSV file
colnames(grouping_mx_complete) <- c("Group No.", 
                               "Group Name", 
                               "Start No.", 
                               "Stop No.",
                               "Start BP",
                               "Stop BP", 
                               "New Members",
                               "Median")
write.csv(grouping_mx_complete[-1, ], file = outfile_grouping_conplete, row.names = FALSE)


