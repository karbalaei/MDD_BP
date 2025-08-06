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


Enrichment_function <-  function(df , name) {
  
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



MDD_sACC = MDD %>% dplyr::filter(FDR.5perc %in% c("sACC" , "Both")) %>% pull(1)

## Define Universe
all_gencode <- c(MDD$geneid , BP$geneid)

head(all_gencode)

all_ensembl <- unique(ss(all_gencode,"\\."))
length(all_ensembl)


all_entrez <- bitr(all_ensembl, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(all_ensembl, fromType = "ENSEMBL", toType = "ENTREZID",  :
#             22.15% of input gene IDs are fail to map...

nrow(all_entrez)
# [1] 14190
u <- all_entrez$ENTREZ


get_signif <- function(  FDR_col,  outFeature, colname = "geneid",  cutoff = 0.05, return_unique = FALSE  ){
    signif <- pull( outFeature[outFeature[ , FDR_col] < cutoff , colname])
    if(return_unique) signif <- unique(signif)
    signif <- signif[!is.na(signif)]
    e <- unique(ss(signif, "\\."))
    entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
    res  <- entrez$ENTREZID
    return(res)
} 

param <- tibble::tribble(
  ~input_vec,         
  "Amygdala", 
  "sACC"
  
)

param <- tibble::tribble(
  ~FDR_col,        ~outFeature, 
  "Amygdala.fdr.p",  MDD, 
  "sACC.fdr.p",  MDD, 
  "Amygdala.fdr.p", BP,
  "sACC.fdr.p",  BP
  )

gene_list = pmap(param ,~get_signif(.x , .y))
#names(gene_list) <-  paste0 (param$FDR_col %>% ss("\\.") , "_" ,  param$outFeature)
     
ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont


plan <- tidyr::expand_grid(
  gene = gene_list,
  ontology = ont
)

names_tibble <- tidyr::expand_grid(
  diagnosis = c("MDD", "BP"),
  brain_region = c("Amygdala", "sACC"),
  ontology = ont
)

names(Enrichment_results) <- apply(names_tibble, 1, paste, collapse = "_")

all_plots <-  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

walk2(names(all_plots), all_plots, ~ggsave(filename = here("graphs" , "twas_Enrichment results" , paste0(.x , ".jpg")), plot = .y, 
                                           height = 14, width = 9))



pdf(here("graphs" , "twas_Enrichment results" , "Enrichment results MDD vs. BP.pdf") , height = 14 , width = 9)

map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results,file = here("results", "enrichGO_twas.rda"))




names(Enrichment_results) <- c(paste0(  names(Enrichment_results), "_", rep(ont)))



Amygdala_MDD = get_signif("Amygdala.fdr.p" , outFeature= MDD)
Amygdala_BP = get_signif("Amygdala.fdr.p" , BP)
sACC_MDD = get_signif("sACC.fdr.p" , MDD)
sACC_BP = get_signif("sACC.fdr.p" , BP)

lapply(param , get_signif)


pwalk(.f = get_signif , .l = param)

signif_genes_A<- (Outlist, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "none" ))


get_signif <- function(outFeature_data, colname = "geneid", FDR_col, cutoff = 0.05, return_unique = FALSE) {
    signif_genes <- pull(outFeature_data[outFeature_data[, FDR_col] < cutoff, colname])
  if (return_unique) {
    signif_genes <- unique(signif_genes)
  }
  signif_genes <- signif_genes[!is.na(signif_genes)]
  
    e <- unique(ss(signif_genes, "\\."))
  
  entrez <- tryCatch({
    bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("Error during ENSEMBL to ENTREZID conversion: ", e$message)
    return(NULL)
  })
  
  if (is.null(entrez) || nrow(entrez) == 0) {
    return(numeric(0)) # Return an empty numeric vector if conversion fails or no IDs
  }
  res <- entrez$ENTREZID
  return(res)
}

run_go_enrichment_pipeline <- function(feature_data_list, universe_genes, cutoff = 0.05) {
  
  param <- tibble::tribble(
    ~FDR_col, ~outFeature_name,
    "Amygdala.fdr.p", "MDD",
    "sACC.fdr.p", "MDD",
    "Amygdala.fdr.p", "BP",
    "sACC.fdr.p", "BP"
  )
  
   gene_sets_af <- param %>%
    mutate(
      gene_list = pmap(., function(FDR_col, outFeature_name) {
        current_outFeature_data <- feature_data_list[[outFeature_name]]
        if (is.null(current_outFeature_data)) {
          warning(paste("Data for", outFeature_name, "not found in feature_data_list. Skipping."))
          return(numeric(0))
        }
        get_signif(
          outFeature_data = current_outFeature_data,
          FDR_col = FDR_col,
          cutoff = cutoff,
          return_unique = TRUE # Assuming you want unique genes for enrichment
        )
      })
    ) %>%
    # Add a name for each gene set for easier identification later
    mutate(list_name = paste(outFeature_name, sub(".fdr.p", "", FDR_col), sep = "_")) %>%
    # Select only the gene_list and name it
    pull(gene_list, name = list_name)
  
  # Remove empty gene sets that might result from no significant genes or conversion issues
  gene_sets_af <- gene_sets_af[lengths(gene_sets_af) > 0]
  
  if (length(gene_sets_af) == 0) {
    message("No significant gene sets were generated. Exiting GO enrichment.")
    return(NULL)
  }
  
  # 3. Prepare for GO enrichment
  ont <- c("BP", "CC", "MF", "ALL")
  names(ont) <- ont
  
  plan <- tidyr::expand_grid(
    gene_list_name = names(gene_sets_af), # Use names to refer back to original lists
    ontology = ont
  )
  
  # 4. Run GO enrichment
  Enrichment_results <- pmap(plan, ~ {
    current_gene_list <- gene_sets_af[[.x]] # Get the actual gene list by name
    current_ontology <- .y
    
    if (length(current_gene_list) == 0) {
      warning(paste("Skipping enrichment for", .x, "with ontology", .y, "due to empty gene list."))
      return(NULL)
    }
    
    tryCatch({
      enrichGO(gene = current_gene_list,
               universe = universe_genes,
               OrgDb = org.Hs.eg.db,
               ont = current_ontology,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.2,
               qvalueCutoff = 0.5,
               readable = TRUE)
    }, error = function(e) {
      message(paste("Error running enrichGO for", .x, "with ontology", current_ontology, ":", e$message))
      return(NULL)
    })
  })
  
  # Name the results for easier access
  names(Enrichment_results) <- paste(plan$gene_list_name, plan$ontology, sep = "_")
  
  # Filter out NULL results (from errors or empty gene lists)
  Enrichment_results <- compact(Enrichment_results)
  
  return(Enrichment_results)
}

all_data <- list(MDD = MDD, BP = BP)


final_enrichment_results <- run_go_enrichment_pipeline(
  feature_data_list = all_data,
  universe_genes = u,
  cutoff = 0.05
          
  gene_set_numbers <-  map_int(final_enrichment_results, nrow)
  