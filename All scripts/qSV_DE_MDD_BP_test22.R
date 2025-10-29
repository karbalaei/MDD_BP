library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(tidyverse)
library(sessioninfo)
library(here)

source(here("codes", "run_DE2.R"))
#source(here("codes", "run_DE_original.R"))

##### Load Data ####

## load rse
load(here("data" , "rse_gene.Rdata"), verbose=TRUE)
load(here("data" , "rse_jxn.Rdata"), verbose=TRUE)
load(here("data" , "rse_exon.Rdata"), verbose=TRUE)
load(here("data" , "rse_tx.Rdata"), verbose=TRUE)


pd = colData(rse_gene)


#### Standardize rowData ####
Reduce(intersect, map(list(rse_gene, rse_exon, rse_jxn, rse_tx), ~colnames(rowData(.x))))
# [1] "meanExprs"    "passExprsCut"

## Add common_feature_id, common_gene_symbol, common_gene_ID
rowData(rse_gene)$common_feature_id <- rowData(rse_gene)$gencodeID
rowData(rse_gene)$common_gene_symbol <- rowData(rse_gene)$Symbol
rowData(rse_gene)$common_gene_id <- rowData(rse_gene)$gencodeID

rowData(rse_exon)$common_feature_id <- rowData(rse_exon)$exon_gencodeID
rowData(rse_exon)$common_gene_symbol <- rowData(rse_exon)$Symbol
rowData(rse_exon)$common_gene_id <- rowData(rse_exon)$gencodeID

rowData(rse_jxn)$common_feature_id <- rownames(rse_jxn)
rowData(rse_jxn)$common_gene_symbol <- rowData(rse_jxn)$Symbol
rowData(rse_jxn)$common_gene_id <- rowData(rse_jxn)$gencodeGeneID
rowData(rse_jxn)$gencodeTx <- unlist(lapply(rowData(rse_jxn)$gencodeTx, toString))

rowData(rse_tx)$common_feature_id <- rowData(rse_tx)$transcript_id
rowData(rse_tx)$common_gene_symbol <- rowData(rse_tx)$gene_name
rowData(rse_tx)$common_gene_id <- rowData(rse_tx)$gene_id

Reduce(intersect, map(list(rse_gene, rse_exon, rse_jxn, rse_tx), ~colnames(rowData(.x))))
# # [1] "meanExprs"          "passExprsCut"       "common_feature_id"  "common_gene_symbol" "common_gene_id"



## qSV data
load(here("data" ,"qSV_mat.Rdata"), verbose = TRUE)



#### Define models ####
regions <- list(sacc = "sACC", amyg = "Amygdala")
region_samples <- map(regions, ~rownames(pd)[pd$BrainRegion == .x])

modSep <- map(regions, ~cbind(model.matrix(~PrimaryDx + AgeDeath + Sex  + mitoRate + 
                                           rRNA_rate +totalAssignedGene + RIN + abs(ERCCsumLogErr) +
                                           snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                                           data= pd[pd$BrainRegion == .x,]),
                              qSV_mat[pd$BrainRegion == .x,]))

map(modSep, colnames)


#### Gene ####
message("\nGENE")
rse_gene_sep <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])

map(rse_gene_sep, dim)


BP_MDD_Gene_Final<- map2(rse_gene_sep, modSep, ~run_DE(rse = .x, model = .y, coef = c("PrimaryDxControl","PrimaryDxBipolar") , save_eBayes = TRUE))

save(BP_MDD_Gene_Final, file = here("results" , "BP_MDD_Gene_qSV_test22.Rdata"), compress=TRUE)

BP_MDD_Gene_Final = transpose(BP_MDD_Gene_Final)

BP_MDD_Gene_Final= BP_MDD_Gene_Final$topTable

#### Exon ####
message("\nEXON")

rse_sep_exon <- map(regions, ~rse_exon[,rse_exon$BrainRegion == .x])
map(rse_sep_exon, dim)

BP_MDD_Exon_Final<- map2(rse_sep_exon, modSep, ~run_DE(rse = .x, model = .y, coef = c("PrimaryDxControl","PrimaryDxBipolar") , save_eBayes = TRUE))
save(BP_MDD_Exon_Final, compress = TRUE, file = here("results" , "BP_MDD_Exon_qSV_test22.Rdata"))

BP_MDD_Exon_Final = transpose(BP_MDD_Exon_Final)

BP_MDD_Exon_Final= BP_MDD_Exon_Final$topTable

#### Jxn ####
message("\nJXN")

rse_sep_jxn <- map(regions, ~rse_jxn[,rse_jxn$BrainRegion == .x])

map(rse_sep_jxn, dim)

BP_MDD_Junction_Final<- map2(rse_sep_jxn, modSep, ~run_DE(rse = .x, model = .y, coef = c("PrimaryDxControl","PrimaryDxBipolar") ,  save_eBayes = TRUE))
save(BP_MDD_Junction_Final, file = here("results" ,"BP_MDD_Junction_qSV_test22.Rdata"))


BP_MDD_Junction_Final = transpose(BP_MDD_Junction_Final)

BP_MDD_Junction_Final= BP_MDD_Junction_Final$topTable


#### Tx ####
message("\nTX")
#load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)
rse_sep_tx <- map(regions, ~rse_tx[,rse_tx$BrainRegion == .x])
map(rse_sep_tx, dim)


BP_MDD_Transcript_Final<- map2(rse_sep_tx, modSep, ~run_DE(.x, .y, run_voom = FALSE, coef = c("PrimaryDxControl","PrimaryDxBipolar") , save_eBayes = TRUE))
save(BP_MDD_Transcript_Final, compress=TRUE, file = here("results" ,"BP_MDD_Transcript_qSV_test22.Rdata"))

rd_tx <- rowData(rse_tx)[,c("common_feature_id","common_gene_symbol","common_gene_id")]

BP_MDD_Transcript_Final = transpose(BP_MDD_Transcript_Final)

BP_MDD_Transcript_Final= BP_MDD_Transcript_Final$topTable

BP_MDD_Transcript_Final<- map_depth(BP_MDD_Transcript_Final, 1, ~cbind(.x, rd_tx[rownames(.x),]))


Outlist_test22 <-  list("gene" = BP_MDD_Gene_Final, "exon" = BP_MDD_Exon_Final, "jxn" = BP_MDD_Junction_Final ,"tx" = BP_MDD_Transcript_Final)


save(Outlist_test22 , file = here("results", "BP_MDD_All_results_test22.RDS"))




## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
