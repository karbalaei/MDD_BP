library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(here)
library(enrichplot)
library(ggnewscale)
library(sessioninfo)
library("KEGGREST")
library(vctrs)
library(tibble)
library(stringr)
library(ggplot2)
library(readxl)

### TWAS ####
MDD <-  read_excel(here("twas/analysis/tables/MDD_BP_Amyg_sACC_FinalOutputTable.xlsx") , sheet = "TWAS Z Scatterplot with FDR and")


BP <-  read_excel(here("twas/analysis_BP/tables/MDD_BP_Amyg_sACC_FinalOutputTable_BP.xlsx") , sheet = "TWAS Z Scatterplot with FDR and")


get_Entrez_id <- function( input_vec){
    input_vec <- unique(input_vec)
    input_vec <- input_vec[!is.na(input_vec)]
    e <- unique(ss(input_vec, "\\."))
    entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
    res  <- entrez$ENTREZ
    return(res)
  } 

#df = MDD


ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont


twas_Enrichment_function <-  function(df , name) {
  
  u <-  get_Entrez_id(df$geneid)

  sACC = df %>% dplyr::filter(FDR.5perc %in% c("sACC" , "Both")) %>% pull(1) 
  #%>% get_Entrez_id()
  Amygdala = df %>% dplyr::filter(FDR.5perc %in% c("Amygdala" , "Both")) %>% pull(1) 
  #%>% get_Entrez_id()
 
  genes <-  list(Amygdala , sACC)
  
  gene_list = map(genes  , get_Entrez_id)
  
  
  plan <- tidyr::expand_grid(
    gene = gene_list,
    ontology = ont
  )
  
  Enrichment_results <-  pmap(plan ,  ~enrichGO(.x,
                                                universe = u, OrgDb = org.Hs.eg.db,
                                                ont = .y , pAdjustMethod = "BH",
                                                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                                readable= TRUE)   )
  
  
  names_tibble <- tidyr::expand_grid(
    brain_region = c("Amygdala", "sACC"),
    ontology = ont
  )
  names(Enrichment_results) <- apply(names_tibble, 1, paste, collapse = "_")
  
  all_plots <-  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))
  
  walk2(names(all_plots), all_plots, ~ggsave(filename = here("graphs" , "twas_Enrichment results" , paste0(.x , "_" , name ,".jpg")), plot = .y, 
                                             height = 14, width = 9))
  
  
  
  pdf(here("graphs" , "twas_Enrichment results" ,paste0("Enrichment results MDD vs. BP_", name,  "_.pdf")) , height = 14 , width = 9)
  
  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))
  
  dev.off()
  
  save(Enrichment_results,file = here("results", paste0(name , "_enrichGO_twas.rda")))
  
  
}

twas_Enrichment_function(MDD , "MDD")
twas_Enrichment_function(BP , "BP")


##### IsoTWAS ####

Amygdala_BP <-  read_excel(here("isotwas/Results/Amygdala_Final_results_BP.xlsx") , sheet = "FDR_filtered_results")

Amygdala_MDD <-  read_excel(here("isotwas/Results/Amygdala_Final_results_MDD.xlsx") , sheet = "FDR_filtered_results")

sACC_BP <-  read_excel(here("isotwas/Results/sACC_Final_results_BP.xlsx") , sheet = "FDR_filtered_results")

sACC_MDD <-  read_excel(here("isotwas/Results/sACC_Final_results_MDD.xlsx") , sheet = "FDR_filtered_results")


u_BP <-  c(read_excel(here("isotwas/Results/Amygdala_Final_results_BP.xlsx") , sheet = "All_results") %>% pull(1) , 
                 read_excel(here("isotwas/Results/sACC_Final_results_BP.xlsx") , sheet = "All_results") %>% pull(1)) %>% unique() %>% get_Entrez_id()

u_MDD <-  c(read_excel(here("isotwas/Results/Amygdala_Final_results_MDD.xlsx") , sheet = "All_results") %>% pull(1) , 
                 read_excel(here("isotwas/Results/sACC_Final_results_MDD.xlsx") , sheet = "All_results") %>% pull(1)) %>% unique()  %>% get_Entrez_id()



isotwas_Enrichment_function <-  function(Amygdala_isotwas , sACC_isotwas,uni ,  name) {
  


  Amygdala_list <-  unique(Amygdala_isotwas$Gene )

    sACC_list <-  unique(sACC_isotwas$Gene)

 
  genes <-  list(Amygdala_list , sACC_list)
  
  gene_list = map(genes  , get_Entrez_id)
  
  
  plan <- tidyr::expand_grid(
    gene = gene_list,
    ontology = ont
  )
  
  Enrichment_results <-  pmap(plan ,  ~enrichGO(.x,
                                                universe = uni, OrgDb = org.Hs.eg.db,
                                                ont = .y , pAdjustMethod = "BH",
                                                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                                readable= TRUE)   )
  
  
  names_tibble <- tidyr::expand_grid(
    brain_region = c("Amygdala", "sACC"),
    ontology = ont
  )
  names(Enrichment_results) <- apply(names_tibble, 1, paste, collapse = "_")
  
  all_plots <-  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))
  
  walk2(names(all_plots), all_plots, ~ggsave(filename = here("graphs" , "isotwas_Enrichment results" , paste0(.x , "_" , name ,".jpg")), plot = .y, 
                                             height = 14, width = 9))
  
  
  
  pdf(here("graphs" , "isotwas_Enrichment results" ,paste0("Enrichment results MDD vs. BP_", name,  "_.pdf")) , height = 14 , width = 9)
  
  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))
  
  dev.off()
  
  save(Enrichment_results,file = here("results", paste0(name , "_enrichGO_isotwas.rda")))
  
  
}



isotwas_Enrichment_function(Amygdala_BP , sACC_BP , u_BP , "BP")

isotwas_Enrichment_function(Amygdala_MDD , sACC_MDD , u_MDD , "MDD")
 

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

  