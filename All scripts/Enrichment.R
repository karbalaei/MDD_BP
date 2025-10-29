### Library #####
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




load(file = here("results" ,"BP_MDD_All_results.RDS"))

## Define Universe
all_gencode <- map_depth(Outlist, 2, "common_gene_id")
map_depth(all_gencode, 2, length)
head(all_gencode$gene$amyg)


length(unlist(all_gencode))
# [1] 1620058

## all ENSEMBL
all_ensembl <- unique(ss(unlist(all_gencode),"\\."))
length(all_ensembl)
# [1] 34635

all_entrez <- bitr(all_ensembl, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
# 31.7% of input gene IDs are fail to map...
nrow(all_entrez)
# [1] 26327
u <- all_entrez$ENTREZ

#### Functions for extracting gene sets ####

get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.05, return_unique = FALSE , UpDown ){
  if (UpDown == "none") {
  signif <- outFeature[[colname]][outFeature$q_PrimaryDxBipolar < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
  } else if ( UpDown == "Up") {
    
    signif <- outFeature[colname][c(outFeature$q_PrimaryDxBipolar < cutoff & outFeature$PrimaryDxBipolar> 0 ), ]

    if(return_unique) signif <- unique(signif)
    signif <- signif[!is.na(signif)]
    return(signif)
    
  } else if (UpDown == "Down") {
    
    
    signif <- outFeature[colname][c(outFeature$q_PrimaryDxBipolar < cutoff & outFeature$PrimaryDxBipolar< 0 ), ]
    if(return_unique) signif <- unique(signif)
    signif <- signif[!is.na(signif)]
    return(signif)
    
  
  }
}




my_flatten <- function (x, use.names = TRUE, classes = "ANY") {
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

my_get_entrez <- function(g){
  e <- unique(ss(g, "\\."))
  entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
  return(entrez$ENTREZ)
}



my_get_entrez_all <- function(g){
  e <- unique(ss(g, "\\."))
  entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
  return(entrez)
}

#### Get signif genes####
signif_genes_whole_list<- map_depth(Outlist, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "none" ))

str(signif_genes_whole_list)

signif_genes_whole_list$All_assigned_genes$amyg =  map(signif_genes_whole_list , 1) %>%  list_c() %>% unique()
signif_genes_whole_list$All_assigned_genes$sacc =  map(signif_genes_whole_list , 2) %>%  list_c() %>% unique()


str(signif_genes_whole_list)



signif_genes_up <- map_depth(Outlist, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "Up" ))

str(signif_genes_up)

signif_genes_down<- map_depth(Outlist, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "Down" ))

str(signif_genes_down)


### Whole  ###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_af <- transpose(signif_genes_whole_list)


str(signif_genes_af )



signif_genes_af_flat <- my_flatten(map_depth(signif_genes_af, 2, unlist))

map_int(signif_genes_af_flat, length)


## Extract unique and convert to entrez
gene_sets_af <- map(signif_genes_af_flat, my_get_entrez)

map_int(gene_sets_af, length)

t(rbind(map_dfr(signif_genes_af_flat, length), map_dfr(gene_sets_af, length)))

#### Run Enrichment ####

ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont


plan <- tidyr::expand_grid(
  gene_list = gene_sets_af,
  ontology = ont
)

Enrichment_results <-  pmap(plan ,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results) <- c(paste0(names(Enrichment_results), "_", rep(ont)))

gene_set_numbers <-  map_int(Enrichment_results, nrow)


to_remove<-  which(gene_set_numbers == 0)

if ( length(to_remove) > 0 ) {

	Enrichment_results <-  Enrichment_results[-(c(to_remove))]

 } 

gene_set_numbers <-  map_int(Enrichment_results, nrow)


gene_set_numbers



all_plots <-  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots), all_plots, ~ggsave(filename = here("graphs" , "Enrichment results" , "whole_list" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))



pdf(here("graphs"  , "Enrichment results" , "whole_list" , "Enrichment results MDD vs. BP.pdf") , height = 14 , width = 9)

map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results,file = here("results", "enrichGO_whole_list.rda"))

### Up ###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_uf <- transpose(signif_genes_up)


str(signif_genes_uf)
# List of 2
# $ Amygdala:List of 4
# ..$ gene: chr [1:86] "ENSG00000189266.11" "ENSG00000142686.7" "ENSG00000117069.14" "ENSG00000179902.12" ...
# ..$ exon: chr [1:1245] "ENSG00000074800.13" "ENSG00000120948.15" "ENSG00000116688.16" "ENSG00000142621.19" ...
# ..$ jxn : chr [1:675] "ENSG00000054523.17" "ENSG00000142621.19" "ENSG00000157191.19" "ENSG00000070831.15" ...
# ..$ tx  : chr [1:5705] "ENSG00000225972.1" "ENSG00000224969.1" "ENSG00000188157.13" "ENSG00000131591.17" ...
# $ sACC    :List of 4
# ..$ gene: chr [1:307] "ENSG00000248527.1" "ENSG00000233954.6" "ENSG00000055070.16" "ENSG00000011009.10" ...
# ..$ exon: chr [1:1981] "ENSG00000248333.8" "ENSG00000008130.15" "ENSG00000157881.13" "ENSG00000116198.12" ...
# ..$ jxn : chr [1:1031] "ENSG00000054523.17" "ENSG00000142657.20" "ENSG00000142655.12" "ENSG00000048707.13" ...
# ..$ tx  : chr [1:1556] "ENSG00000237094.11" "ENSG00000248527.1" "ENSG00000131591.17" "ENSG00000175756.13" ...


signif_genes_uf_flat <- my_flatten(map_depth(signif_genes_uf, 2, unlist))
map_int(signif_genes_uf_flat, length)


# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene
# 86          1245           675          5705           307
# sACC.exon      sACC.jxn       sACC.tx
# 1981          1031          1556


## Extract unique and convert to entrez
gene_sets_uf <- map(signif_genes_uf_flat, my_get_entrez)

# 1: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   12.79% of input gene IDs are fail to map...
# 2: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   4.18% of input gene IDs are fail to map...
# 3: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   1.93% of input gene IDs are fail to map...
# 4: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   14.07% of input gene IDs are fail to map...
# 5: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   7.49% of input gene IDs are fail to map...
# 6: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   4.14% of input gene IDs are fail to map...
# 7: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   3.01% of input gene IDs are fail to map...
# 8: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   14.91% of input gene IDs are fail to map...


map_int(gene_sets_uf, length)

# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene
# 76          1216           665          5124           287
# sACC.exon      sACC.jxn       sACC.tx
# 1928          1005          1384



t(rbind(map_dfr(signif_genes_uf_flat, length), map_dfr(gene_sets_uf, length)))

#                [,1] [,2]
# Amygdala.gene   86   76
# Amygdala.exon 1245 1216
# Amygdala.jxn   675  665
# Amygdala.tx   5705 5124
# sACC.gene      307  287
# sACC.exon     1981 1928
# sACC.jxn      1031 1005
# sACC.tx       1556 1384


#### Run Enrichment ####

plan_up <- tidyr::expand_grid(
  gene_list = gene_sets_uf,
  ontology = ont
)

Enrichment_results_up <-  pmap(plan_up,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results_up) <- c(paste0(names(Enrichment_results_up), "_", rep(ont)))

class(Enrichment_results_up)


gene_set_numbers_up <-  map_int(Enrichment_results_up, nrow)


gene_set_numbers_up



to_remove_up<-  which(gene_set_numbers_up == 0)

if ( length(to_remove_up) > 0 ) {

	Enrichment_results_up <-  Enrichment_results_up[-(c(to_remove_up))]

 } 

gene_set_numbers_up <-  map_int(Enrichment_results_up, nrow)


gene_set_numbers_up




all_plots_up <-  map2(Enrichment_results_up, names(Enrichment_results_up), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots_up), all_plots_up, ~ggsave(filename = here("graphs" , "Enrichment results" , "up" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))




pdf(here("graphs" , "Enrichment results" , "up" , "Enrichment results MDD vs. BP_up.pdf") , height = 14 , width = 9)

map2(Enrichment_results_up, names(Enrichment_results_up), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()


save(Enrichment_results_up,file = here("results", "enrichGO_up.rda"))


### down###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_df <- transpose(signif_genes_down)


str(signif_genes_df)

# $ Amygdala:List of 4
# ..$ gene: chr [1:249] "ENSG00000078369.17" "ENSG00000162585.16" "ENSG00000097021.19" "ENSG00000132906.17" ...
# ..$ exon: chr [1:1847] "ENSG00000228794.8" "ENSG00000127054.19" "ENSG00000221978.11" "ENSG00000197530.12" ...
# ..$ jxn : chr [1:782] "ENSG00000160072.19" "ENSG00000197530.12" "ENSG00000158109.14" "ENSG00000130939.18" ...
# ..$ tx  : chr [1:395] "ENSG00000278791.1" "ENSG00000264293.2" "ENSG00000078369.17" "ENSG00000236963.6" ...
# $ sACC    :List of 4
# ..$ gene: chr [1:239] "ENSG00000131584.18" "ENSG00000127054.19" "ENSG00000197530.12" "ENSG00000078369.17" ...
# ..$ exon: chr [1:1565] "ENSG00000188976.10" "ENSG00000188157.13" "ENSG00000078808.16" "ENSG00000131584.18" ...
# ..$ jxn : chr [1:696] "ENSG00000197530.12" "ENSG00000158109.14" "ENSG00000048707.13" "ENSG00000189337.16" ...
# ..$ tx  : chr [1:551] "ENSG00000278791.1" "ENSG00000078808.16" "ENSG00000131584.18" "ENSG00000264293.2" ...



signif_genes_df_flat <- my_flatten(map_depth(signif_genes_df, 2, unlist))
map_int(signif_genes_df_flat, length)
 
# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene
# 249          1847           782           395           239
# sACC.exon      sACC.jxn       sACC.tx
# 1565           696           551


## Extract unique and convert to entrez
gene_sets_df <- map(signif_genes_df_flat, my_get_entrez)


# 1: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   8.03% of input gene IDs are fail to map...
# 2: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   4.11% of input gene IDs are fail to map...
# 3: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   1.92% of input gene IDs are fail to map...
# 4: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   16.46% of input gene IDs are fail to map...
# 5: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   7.53% of input gene IDs are fail to map...
# 6: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   3.83% of input gene IDs are fail to map...
# 7: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.01% of input gene IDs are fail to map...
# 8: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   14.52% of input gene IDs are fail to map...



map_int(gene_sets_df, length)

# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene
# 441          1995           772          2886           433
# sACC.exon      sACC.jxn       sACC.tx
# 1936           685          3028

t(rbind(map_dfr(signif_genes_df_flat, length), map_dfr(gene_sets_df, length)))

#                [,1] [,2]
# Amygdala.gene  249  441
# Amygdala.exon 1847 1995
# Amygdala.jxn   782  772
# Amygdala.tx    395 2886
# sACC.gene      239  433
# sACC.exon     1565 1936
# sACC.jxn       696  685
# sACC.tx        551 3028

#### Run Enrichment ####


plan_down <- tidyr::expand_grid(
  gene_list = gene_sets_df,
  ontology = ont
)

Enrichment_results_down<-  pmap(plan_down ,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results_down) <- c(paste0(names(Enrichment_results_down), "_", rep(ont)))

class(Enrichment_results_down)


gene_set_numbers_down<-  map_int(Enrichment_results_down, nrow)


gene_set_numbers_down

to_remove_down<-  which(gene_set_numbers_down == 0)

if ( length(to_remove_down) > 0 ) {

	Enrichment_results_down <-  Enrichment_results_down[-(c(to_remove_down))]

 } 

gene_set_numbers_down <-  map_int(Enrichment_results_down, nrow)


gene_set_numbers_down

all_plots_down <-  map2(Enrichment_results_down, names(Enrichment_results_down), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots_down), all_plots_down, ~ggsave(filename = here("graphs" , "Enrichment results" , "down" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))




pdf(here("graphs" , "Enrichment results" , "down"  , "Enrichment results MDD vs. BP_down.pdf") , height = 14 , width = 9)

map2(Enrichment_results_down, names(Enrichment_results_down), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()


save(Enrichment_results_down,file = here("results", "enrichGO_down.rda"))






## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

