#!/bin/bash
## These mkdir steps + ln -s + "mkdir -p logs/Amygdala_gene" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## To avoid having to change the file permissions later
## From https://twitter.com/fellgernon/status/1258455434073124865?s=20
## and https://www.cyberciti.biz/tips/understanding-linux-unix-umask-value-usage.html
umask u=rwx,g=rwx,o= ## equivalent to umask 007

## Required order for running this code:
cd ../twas/Amygdala_gene

find ./expression_prediction -type f -name "*.profile" -exec basename {} \; > ./Amygdala_gene_list_prediced.txt.profile


library(data.table)

# --- 1. Define Paths and Files ---
file_list_path <- "Amygdala_gene_list_prediced.txt.profile"
output_file <- "merged_prediced_expression_Amygdala.tsv"

INPUT_DIR <- "./expression_prediction/"

# Read the list of files to merge
files_to_merge <- readLines(file_list_path)

# --- 2. Function to Extract Gene Name ---
# Extracts "gene_9710" from "gene_9710_prediction.profile"
extract_gene_name <- function(filename) {
  # Remove everything from "_prediction.profile" to the end
  sub("_prediction\\.profile$", "", filename)
}

# --- 3. Iterate and Merge ---

# Read the first file to initialize the master data table
cat("Reading and initializing with:", files_to_merge[1], "\n")
master_dt <- fread(files_to_merge[1], select = c(1, 2, 6)) # Select Col 1, 2, and 6
setnames(master_dt, c("V1", "V2", "V6"), c("ID1", "ID2", extract_gene_name(files_to_merge[1])))

# Iterate through the rest of the files (starting from the second file)
for (i in 2:length(files_to_merge)) {
  current_file_base <- files_to_merge[i]
  current_file_path <- paste0(INPUT_DIR, current_file_base) # CONSTRUCT FULL PATH
  current_gene <- extract_gene_name(current_file_base)
  cat("Merging:", current_file_path, "\n")
  
  # Read the current file, selecting only the necessary columns
  current_dt <- fread(current_file_path, select = c(1, 2, 6))
  # Rename columns immediately for merging
  setnames(current_dt, c("ID1", "ID2", current_gene))
  
  # Perform the merge
  master_dt <- merge(master_dt, current_dt, 
                     by = c("ID1", "ID2"), all = TRUE)
}
# --- 4. Write Final Output ---
fwrite(master_dt, file = output_file, sep = "\t", na = "NA")

cat("Merge complete. Output saved to:", output_file, "\n")

# Define the output directory
OUTPUT_DIR="./expression_prediction"
mkdir -p ${OUTPUT_DIR}

#Rscript make_expression_prediction.R [wgt.RDat file] > [expression_prediction_FILE] 

#plink --bfile [GENOTYPES] --expression_prediction [expression_prediction_FILE] 1 2 4


split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.profile ./Amygdala_gene_list_expression_prediction.profile ./Amygdala_gene_list_expression_prediction_

cd ../sACC_gene

find ./expression_prediction -type f -name "*.wgt.RDat" -exec basename {} \; > ./sACC_gene_list_prediced.txt


#Rscript make_expression_prediction.R [wgt.RDat file] > [expression_prediction_FILE] 

#plink --bfile [GENOTYPES] --expression_prediction [expression_prediction_FILE] 1 2 4


split -l 5000 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.profile ./sACC_gene_list_expression_prediction.profile ./sACC_gene_list_expression_prediction_

# Define the output directory
OUTPUT_DIR="./expression_prediction"
mkdir -p ${OUTPUT_DIR}
