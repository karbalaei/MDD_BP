
library("SummarizedExperiment")
library("here")
library(purrr)
library(tidyverse)

options(scipen=999)


load(here("data" ,'rse_jxn.Rdata'), verbose=TRUE)


ind <-  which(rse_jxn$PrimaryDx=="MDD" | rse_jxn$PrimaryDx=="Bipolar")

rse_jxn <-  rse_jxn[ ,ind]  

rse_jxn$PrimaryDx <- droplevels(rse_jxn$PrimaryDx)
rse_jxn$PrimaryDx <- relevel(rse_jxn$PrimaryDx, "Bipolar")



junc_files <- gsub("_accepted_hits.sorted.bam","_junctions_primaryOnly_regtools.bed", gsub("HISAT2_out","Counts/junction",
                    gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_jxn$bamFile)))

junc_files <- gsub("hydeZandi_LIBD2005/goesHyde_mdd_rnaseq","hydeGoes_LIBD3010/goesHyde_mdd_rnaseq", junc_files)


leafcutter_ids = gsub("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/preprocessed_data/Counts/junction/", "" ,   junc_files)
leafcutter_ids = gsub("/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/zandiHyde_bipolar_rnaseq/preprocessed_data/Counts/junction/", "" ,   leafcutter_ids)

rse_jxn$leafcutter_ids = leafcutter_ids

pd = colData(rse_jxn) %>% as.data.frame() %>%  rownames_to_column("id")


load(here("data" ,"qSV_mat.Rdata"), verbose = TRUE)


qSV_mat = qSV_mat[pd$id ,]


regions <- list(sacc = "sACC", amyg = "Amygdala")


modSep <- map(regions, ~cbind(model.matrix(~PrimaryDx + AgeDeath + Sex  + mitoRate + 
                                           rRNA_rate +totalAssignedGene + RIN + abs(ERCCsumLogErr) +
                                           snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                                           data= pd[pd$BrainRegion == .x,]),
                              qSV_mat[pd$BrainRegion == .x,] , id = pd[pd$BrainRegion == .x , "leafcutter_ids"]))

map(modSep, colnames)

write.table(modSep$sacc, file= here("data" , "modSep_sACC_MDDBP.txt") , row.names = F , sep = "\t" , , quote=F)
write.table(modSep$amyg, file= here("data" , "modSep_Amygdala_MDDBP.txt") , row.names = T , sep = "\t" ,  quote=F)



Amygdala_model = read.delim(file= here("data" , "modSep_Amygdala_MDDBP.txt") , sep = "\t")
sACC_model = read.delim(file= here("data" , "modSep_sACC_MDDBP.txt") , sep = "\t")


sACC_model$PrimaryDxMDD = if_else(sACC_model$PrimaryDxMDD =="1" , "MDD" , "BP")
sACC_model$SexM = if_else(sACC_model$SexM =="0" , "Female" , "Male")


Amygdala_model$PrimaryDxMDD = if_else(Amygdala_model$PrimaryDxMDD =="1" , "MDD" , "BP")
Amygdala_model$SexM = if_else(Amygdala_model$SexM =="0" , "Female" , "Male")

Amygdala_model[ , -1] %>% relocate(id) %>% 
write.table(file= here("data" , "modSep_Amygdala_for_leafcutter.txt") , row.names = F , sep = "\t" ,  quote=F , col.names = F)

sACC_model[ , -1] %>% relocate(id) %>% 
write.table( file= here("data" , "modSep_sACC_for_leafcutter.txt") , row.names = F , sep = "\t" ,  quote=F ,  col.names = F)



## Save Models



junc_files <- gsub("_accepted_hits.sorted.bam","_junctions_primaryOnly_regtools.bed", gsub("HISAT2_out","Counts/junction",
                    gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_jxn$bamFile)))

junc_files <- gsub("hydeZandi_LIBD2005/goesHyde_mdd_rnaseq","hydeGoes_LIBD3010/goesHyde_mdd_rnaseq", junc_files)


leafcutter_ids = gsub("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/preprocessed_data/Counts/junction/", "" ,   junc_files)
leafcutter_ids = gsub("/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/zandiHyde_bipolar_rnaseq/preprocessed_data/Counts/junction/", "" ,   leafcutter_ids)


# identify the folders
new.folder <- "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/leafcutter/junc_files"

# find the files that you want

# copy the files to the new folder
file.copy(junc_files, new.folder)


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()




