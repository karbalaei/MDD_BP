## Load plink before starting R:
# module load plink/1.90b6.6
## Also load twas fusion code
# module load fusion_twas/github

library("SummarizedExperiment")
library("jaffelab")
library("data.table")
library("sessioninfo")
library("getopt")
library("BiocParallel")
library("tidyr")
library("here")
library(readr)
library(dplyr)
library(stringr)
library(tibble)

## For styling this script
# library("styler")
# styler::style_file("build_bims.R", transformers = biocthis::bioc_style())

## Without this, the memory use blows up
## getDTthreads() will detect 64 threads in some cases here
setDTthreads(threads = 1)

## Flags that are supplied with RScript
## Specify parameters
## Tissues: Amygdala, SACC
spec <- matrix(
  c(
    'region',
    'r',
    1,
    'character',
    'Either Amygdala or SACC',
    'cores',
    'c',
    1,
    'integer',
    'Number of cores to use. Use a small number',
    'help' ,
    'h',
    0,
    'logical',
    'Display help'
  ),
  byrow = TRUE,
  ncol = 5
)
opt <- getopt(spec)

## For an interactive test
if (FALSE) {
  opt <-
    list(
      'region' = 'Amygdala',
      "feature" = "transcript",
      "cores" = 1,
      "test" = TRUE
    )
}
#opt$region <- "Amygdala"
# opt$cores <- "cores"


# opt$region <- tolower(opt$region)
opt$feature <- "transcript"

stopifnot(opt$region == "Amygdala" | opt$region == "sACC")

dir.create(paste0(opt$region, "_", opt$feature), showWarnings = FALSE)

# default arguments for flags
if (is.null(opt$test)) {
  opt$test <- FALSE
}
if (is.null(opt$cores)) {
  ## Default to 1 core
  opt$cores <- 1
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

## Show the options used
message(paste(Sys.time(), "options used"))
print(opt)

# feat <- opt$feature
# reg <- opt$region

# official
load_rse <- function(feat, reg) {
  message(paste(Sys.time(), "loading expression data"))
  
  # expmnt data
  load(Sys.glob(
    here::here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , 
               paste0(opt$region, "_rda"),
               paste0(opt$region, "_hg38_rseTranscript_rawCounts_allSamples_n*.Rdata")
    )
  ), verbose = TRUE)
  ## Could be more complicated later on for exon, jxn, tx
  rse <- rse_tx
  assays(rse)$raw_expr <- assays(rse_tx)$tpm
  
  load(here::here(
    "hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , 
    paste0(opt$region, "_rda"),
    "TranscriptPCs.Rdata"
  ),
  verbose = TRUE)
  pcs <- TranscriptPCs
  
  pd = colData(rse)
  
  colnames(pd)
  
  colData(rse) <- cbind(colData(rse), pcs)
  print(head(rse))
  
  mod <-
    model.matrix( ~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + pcs,
                  data = colData(rse))
  
  message(paste(Sys.time(), 'cleaning expression'))
  assays(rse) <- List('raw_expr' = assays(rse)$raw_expr,
                      'clean_expr' = cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 3))
  
  ## Regress out effects. If we had a diagnosis variable (Dx), we would use it
  ## first, then use P = 2 in cleaningY()
  
  message(paste(Sys.time(), "switch column names to BrNum"))
  stopifnot(!any(duplicated(colnames(rse))))
  colnames(rse) <- colData(rse)$BrNum
  
  return(rse)
  
}

dir.create(paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , opt$region, "_", opt$feature), showWarnings = FALSE)


rse_file <-
  file.path("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas", paste0(opt$region, "_", opt$feature), 'working_rse.RData')

if (file.exists(rse_file) == FALSE) {
  rse <- load_rse(here(
    "hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" ,opt$feature, opt$region))
  message(paste(Sys.time(), "saving the rse file for later at", rse_file))
  save(rse, file = rse_file)
} else {
  message(paste(Sys.time(), "loading previous rse file", rse_file))
  load(rse_file, verbose = TRUE)
}
message(Sys.time(), " working RSE dimensions")
rse_tx = rse
print(dim(rse_tx))

rm(rse)
heritability_info = read.table(file = here("hydeGoes_scSeq_mdd/MDD_vs_BP/twas" ,  paste0("All_hsq_" , opt$region, ".tsv")) , sep=" ") %>% dplyr::select(-5) #from twas analysis

colnames(heritability_info) = c("File_name" , "Heritability" , "LRT" , "pvalue")

heritability_filtered = heritability_info %>% dplyr::filter(Heritability > 0  & pvalue <= 0.05) 


Gene_index = str_remove_all(heritability_filtered$File_name  , "out_files/gene_") %>% as.numeric() %>% sort()

Gene_address =  paste0("./Amygdala_gene/bim_files/Amygdala_gene_" , Gene_index) 

write.table(Gene_index , file= here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , paste0(opt$region, "_", opt$feature ,"/selected_gene_index_" , opt$region , ".txt")) , row.names= F , col.names = F)

write.table(Gene_address , file= here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , paste0(opt$region, "_", opt$feature ,"/selected_gene_address_" , opt$region , ".txt")), row.names= F , col.names = F , quote=F)

rse_file_gene <-
  file.path("hydeGoes_scSeq_mdd/MDD_vs_BP/twas", paste0(opt$region, "_gene"), 'subsetted_rse.RData')

message(paste(Sys.time(), "loading previous rse gene file from twas", rse_file_gene))
load(rse_file_gene, verbose = TRUE)

rse_gene = rse
message(Sys.time(), " working RSE dimensions")
print(dim(rse_gene))

rm(rse)
rowData(rse_gene)$common_feature_id <- rowData(rse_gene)$gencodeID
rowData(rse_gene)$common_gene_symbol <- rowData(rse_gene)$Symbol
rowData(rse_gene)$common_gene_id <- rowData(rse_gene)$gencodeID


rowData(rse_tx)$common_feature_id <- rowData(rse_tx)$transcript_id
rowData(rse_tx)$common_gene_symbol <- rowData(rse_tx)$gene_name
rowData(rse_tx)$common_gene_id <- rowData(rse_tx)$gene_id


gene_isoforms_count = rowData(rse_tx) %>% as.data.frame() %>% dplyr::select(c("common_feature_id" , "common_gene_id")) %>% group_by(common_gene_id) %>% dplyr::summarize( n = n_distinct(common_feature_id)) %>% dplyr::filter(n > 1) %>% pull(1)

#gene_list_multi_tx =  rowData(rse_gene) %>% as.data.frame() %>% dplyr::filter(NumTx > 1) %>% dplyr::pull(common_gene_id)

gene_name_index = data.frame(common_gene_id =  rowData(rse_gene)$common_gene_id)  %>% rownames_to_column("index") %>% relocate("index")

gene_list_isotwas = gene_name_index %>% dplyr::filter(common_gene_id %in% gene_isoforms_count & index %in% Gene_index)

#gene_list_filtered = gene_name_index %>% dplyr::filter(common_gene_id %in% gene_list_multi_tx & index %in% Gene_index)

gene_details = rowRanges(rse_gene) %>% as.data.frame() %>% dplyr::select(c(seqnames , common_gene_id))

#gene_list_isotwas_details = gene_list_filtered %>% left_join(gene_details , by = "common_gene_id")

gene_list_isotwas_details = gene_list_isotwas %>% left_join(gene_details , by = "common_gene_id")

gene_list_isotwas_details$seqnames =  str_remove_all(gene_list_isotwas_details$seqnames , pattern = "chr") %>% as.numeric()

write.table(gene_list_isotwas_details , file= here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , paste0(opt$region, "_", opt$feature ,"/gene_list_isotwas_" , opt$region , ".txt")) , row.names= F , col.names = T , quote=F , sep = "\t")

rse_tx = rse_tx[rowData(rse_tx)$common_gene_id %in% gene_list_isotwas$common_gene_id , ] 

expression_clean = t(assays(rse_tx)$clean_expr)

rse_info = colData(rse_tx) %>% as.data.frame() %>% rownames_to_column( var = "SampleID") 

stopifnot(identical(rownames(expression_clean), rse_info$SampleID))

rownames(expression_clean) =   str_replace_all( rse_info$genoSample , "_" , ":")

rowData(rse_tx) %>% as.data.frame() %>% dplyr::select(c("common_feature_id" , "common_gene_id")) %>%
  write.table(file= here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , paste0(opt$region, "_", opt$feature ,"/gene_transcript_" , opt$region , ".txt")) , row.names= F , col.names = T , quote=F , sep="\t")


save(expression_clean, file = here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas", paste0(opt$region, "_", opt$feature , "/expression_matrix.rda")))


## Using the files where the SNV names have been made unique
## using make.names(unique = TRUE)
## Details at filter_data/filter_snps.R


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
