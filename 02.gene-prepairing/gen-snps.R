# R program for generate Case-Control SNPs from real SNPs data
# Code by Sirikanlaya Sookkhee

# Edited 27 AUG 2016 22:03 PM
# By Sirikanlaya Sookkhee

# Edited 30 AUG 2016 06:30 PM
# By Fazil Baksh, Pianpool Kirdwichai and Sirikalaya Sookkhee 

# *****************************************************************************#

# Funciton for generate only one new genotype
# @param genotype : matrix for store genotype data from CSV file
# @param dis_snp_pos : position of SNPs for predict 'Case' or 'Control'
# @param loop_count : use for set seed number
# @return new_genotype
generateGenotype <- function(genotypes, dis_snp_pos, loop_count) {
     # Generate new genotype from real data
     # @gen: Step 1 - random 2 samples form genotypes
     set.seed(as.numeric(Sys.time()) + loop_count)
     
     sample_total <- nrow(genotypes)
     
     # random individual's index
     samples_index <- sample(1:sample_total, 2, replace = T)
     
     # @gen: Step 2 - calculate x_spec
     x_spec <- as.numeric(genotypes[samples_index[1], dis_snp_pos]
                          + genotypes[samples_index[2], dis_snp_pos])
     
     # @gen: Step 3 - Find probabiliy value using logistic regression
     # Edited 2 SEP 2016 6:05 PM
     alpha <- 0.0
     gene_effect <- 0.0
     
     exp_value <- exp(alpha + (gene_effect * x_spec))
     prob <- exp_value / (1 + exp_value)
     
     # @gen: Step 4 - Using uniform random to cutoff \
     #                new genotype is 'Case' or 'Control'
     cutoff_value <- runif(1, 0, 1)
     if (cutoff_value < prob) {
          disease <- 1
     } else {
          disease <- 0
     }
     
     # @gen: Step 5 - create new Individual by add vector A and vector B together
     tmp_genotype <- genotypes[samples_index[1],] + genotypes[samples_index[2],]
     
     # @gen: step 6 - insert disease or Y value into fist index of list
     new_genotype <- append(tmp_genotype, disease, 0)
     new_genotype <- append(new_genotype, samples_index[1], 1)
     new_genotype <- append(new_genotype, samples_index[2], 2)
     
     # @gen : step 7 - return only one new genotype
     return(new_genotype)
}

#--------------------------------------------------------------------
# Main Program

# Start program
# by get command line argument
# arg[1] is start number of replicated
# arg[2] is stop number of replicated
replicate_range <- commandArgs(TRUE)

if (length(replicate_range) == 2 && replicate_range[1] < replicate_range[2]) {
     # @Section 1 : Configure values
     setwd("~/ResearchCode/SNPsR")
     snpsName_path <- "input/gen-snps/NameONLY13479.csv"
     genotype_path <- "input/gen-snps/New_data_3008persons_13479SNPs.csv"
     
     dis_snp = 'rs3789038'
     number_of_population <- 3000
     number_of_case <- number_of_population / 2
     
     # Read input data from CSV files
     snpsName_data <- read.csv(snpsName_path, header = TRUE)
     genotypes_data <- read.csv(genotype_path, header = TRUE)
     
     # Get number of SNPs from genotype data frame
     snps_total <- ncol(genotypes_data)
     
     # Convert snpsName data frame to vector
     # And find SNPs position
     snpsName <- snpsName_data[,1]
     dis_snp_pos <- match(dis_snp, snpsName)
     
     # Convert genotype data frame to matrix
     genotypes <- data.matrix(genotypes_data)
     
     # @Section 2 : Generate SNPs
     for (i in replicate_range[1]:replicate_range[2]) {
          # Create matrix for store 'Case' and 'Control'
          # 1 - 500 : is a Case
          # 501 - 1000 is a Control
          all_data <- matrix(0, nrow = number_of_population, ncol = snps_total + 3)
          case_index  <- 0
          control_index <- number_of_case
          
          loop_count <- 1
          while (case_index < number_of_case || 
                 control_index < number_of_population) {
               
               all_data_index <- -1
               
               # Call function generateGenotype
               # funtion generateGenotype() will return new genotype as a List
               new_genotype <- generateGenotype(genotypes, dis_snp_pos, loop_count)
               loop_count = loop_count + 1
               
               # Classify genotypes between 'Case' and 'Control' by using Y value
               # y_value = 1 is 'Case'
               # y_value = 0 is 'Control'
               y_value <- new_genotype[1]
               
               if (y_value == 1 && case_index < number_of_case) {
                    case_index <- case_index + 1
                    all_data_index <- case_index
               } 
               else if (y_value == 0 && control_index < number_of_population) {
                    control_index <- control_index + 1
                    all_data_index <- control_index
               }
               
               # Store new gentoype/individual with Y value into all_data
               # and classify it is 'Case' or 'Control'
               if (all_data_index != -1) {
                    all_data[all_data_index, ] <- new_genotype
               }
          }
          
          # Write Case & Control data set to CSV file ******************************#
          # Set result's filename by using data and time format
          outfile_path <- "/Volumes/Sirikanlaya/gen-snps-phisan/0.0/"
          outfile_case_control <- paste(outfile_path,
                                        i, 
                                        ".csv", 
                                        sep = "")
          
          # Set header to first row of result matrix
          colnames(all_data) <- c("Y value", 
                                  "Parent 1", 
                                  "Parent 2", 
                                  paste("SNP", 1:snps_total))
          
          # Write all data to CSV
          write.csv(all_data, file = outfile_case_control, row.names = FALSE)
          
          cat(sprintf("Replicated #%s Already created, file %s\n", 
                      i, outfile_case_control))
     }
} else {
     print("Please, Enter replicated range for example 1 100 !!")
}