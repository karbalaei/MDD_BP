library(data.table)
library(here)


# --- 1. Define Paths and Files ---
file_list_path <- here("twas/sACC_gene" , "sACC_gene_list_prediced.txt")
output_file <- here("twas/sACC_gene","merged_prediced_expression_sACC.tsv")

INPUT_DIR <- here("twas/sACC_gene", "expression_prediction")

# Read the list of files to merge
files_to_merge <- readLines(file_list_path)

# --- 2. Function to Extract Gene Name ---
# Extracts "gene_i" from "gene_i_prediction.profile"
extract_gene_name <- function(filename) {
  # Remove everything from "_prediction.profile" to the end
  sub("_prediction\\.profile$", "", filename)
}

# --- 3. Iterate and Merge ---

# Read the first file to initialize the master data table
cat("Reading and initializing with:", files_to_merge[1], "\n")
master_dt <- fread(here("twas/sACC_gene/expression_prediction", files_to_merge[1]), select = c(1, 2, 6)) # Select Col 1, 2, and 6
setnames(master_dt, c("ID1", "ID2", extract_gene_name(files_to_merge[1])))

# Iterate through the rest of the files (starting from the second file)
for (i in 2:length(files_to_merge)) {
  #current_file_base <- files_to_merge[i]
  #current_file_path <- c(INPUT_DIR, current_file_base) # CONSTRUCT FULL PATH
  current_gene <- extract_gene_name(files_to_merge[i])
  #name(current_file_base)
  cat("Merging:", current_gene, "\n")
  
  # Read the current file, selecting only the necessary columns
  current_dt <- fread(here("twas/sACC_gene/expression_prediction", files_to_merge[i]),
					select = c(1, 2, 6))
  # Rename columns immediately for merging
  setnames(current_dt, c("ID1", "ID2", current_gene))
  
  # Perform the merge
  master_dt <- merge(master_dt, current_dt, 
                     by = c("ID1", "ID2"), all = TRUE)
}
# --- 4. Write Final Output ---
fwrite(master_dt, file = output_file, sep = "\t", na = "NA")

cat("Merge complete. Output saved to:", output_file, "\n")


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
