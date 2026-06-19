library(here)
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)

# Load functions from three separate files
source("make_gene_plot_Final2.R")
source("make_cluster_table_Final.R")
source("make_cluster_plot_Final.R")

load(here("leafcutter/covariates_results/step_5_Amygdala.Rdata"))

load(file = here("results" ,"BP_MDD_All_results.RDS"))

DEGs <-   Outlist$gene$amyg %>% dplyr::filter(q_PrimaryDxBipolar< 0.05) %>% pull(common_gene_symbol) %>% unique()

clusters <- clusters %>% dplyr::mutate(gene = gsub("<i>|</i>", "", gene, fixed = FALSE)) %>% 
  dplyr::mutate(DEGs = if_else(gene %in% DEGs , "yes" , "no"))


message("creating tables of clusters")

walk(clusters$clusterID, ~filter_intron_table(introns, .x, toSave=TRUE , region = "Amygdala")) 
 
message("creating graphs of clusters")
 
walk( clusters$clusterID , ~make_cluster_plot(cluster_to_plot = .x ,exons_table = exons_table ,cluster_ids = cluster_ids , introns =  introns , meta = meta , counts = counts ,snp_pos <-  NA , file_path = here("graphs" ,  "leafcutter" , region= "Amygdala" ,  paste0(.x , "_graph.pdf"))))

message("creating graphs of genes")

clusters_gene <- clusters %>% 
  dplyr::filter(gene !=".")

walk2( clusters_gene$gene , clusters_gene$clusterID , ~make_gene_plot(gene_name =  .x ,clusterID <- .y , cluster_list = clusters , introns = introns , introns_to_plot = introns_to_plot , exons_table = exons_table , snp_pos <-  NA , min_exon_length <- 0.5  , region = "Amygdala"))



write.table(clusters, file = here("results" , "leafcutter" , paste0("Amygdala" , "_sig_clusters.csv")) , sep = ",", 
            quote = FALSE, row.names = FALSE)

# Introns
intron_df <- introns %>%
  mutate(coord = paste0(chr, ":", start, "-", end)) %>%
  mutate(ensemblID = gsub("_[0-9]$", "", ensemblID) ) %>%
  dplyr::select(clusterID, gene, ensemblID, coord, verdict, deltaPSI = deltapsi)


write.table(intron_df, file = here("results" , "leafcutter" , paste0("Amygdala" , "_sig_introns.csv")), sep = ",", 
            quote = FALSE, row.names = FALSE)




## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
