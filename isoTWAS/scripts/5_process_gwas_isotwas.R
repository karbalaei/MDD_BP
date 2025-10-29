library(data.table)
library(dplyr)

setDTthreads(1)

# Read in Fernando's GWAS
hg38_gwas <-
    fread(
        "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/PGC_UKB_depression_genome-wide_hg38.txt"
    )

# reorder columns
col_order <- c("chr",
               "MarkerName",
               "bp",
               "A1",
               "A2",
               "Freq",
               "LogOR",
               "StdErrLogOR",
               "P")

hg38_gwas <- hg38_gwas[, ..col_order]

names(hg38_gwas)[1:3] <- c("CHR", "SNP", "BP")

# read bim file with unique rsIDs
uniq_bim <-
    fread(
        "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/twas/filter_snps/Amygdala_unique_snps_bim/mdd_bpd_maf01.rsidAmygdala_uniqueSNPs.bim")
    
colnames(uniq_bim) <- c("CHR", "SNP", "dummy", "BP", "A1", "A2")

# check for overlap between gwas and unique bim file
uniq_bim$new_test <- paste0(uniq_bim$CHR, "_", uniq_bim$BP)

hg38_gwas$new_test <-
    paste0(gsub("chr", "" , hg38_gwas$CHR), "_", hg38_gwas$BP)

table(hg38_gwas$new_test %in% uniq_bim$new_test)

# FALSE    TRUE 
# 823651 7656313 

# hg38_bb <- head(hg38_gwas)
# bim_bb <- head(uniq_bim)

# merge in unique rsIDs
hg38_gwas <-
    merge(hg38_gwas[,-"SNP"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas$new_test <- NULL

col_order <- c(
    "CHR",
    "SNP",
    "BP",
    "A1",
    "A2",
    "Freq",
    "LogOR",
    "StdErrLogOR",
    "P")

hg38_gwas <- hg38_gwas[, ..col_order]

# hg38_gwas$effect <- log(hg38_gwas$OR)

hg38_gwas$Z <- hg38_gwas$LogOR / hg38_gwas$StdErrLogOR


hg38_gwas$N <- 1/((hg38_gwas$StdErrLogOR^2)*(2*hg38_gwas$Freq*(1-hg38_gwas$Freq)))


setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/hg38")

hg38_gwas_isotwas <-
    hg38_gwas[, c("SNP", "A1", "A2", "Z" , "N" , "CHR" , "BP" , "StdErrLogOR")]


colnames(hg38_gwas_isotwas) = c("SNP", "A1", "A2", "Z" , "N" , "Chromosome" , "Position" , "SE")

# colnames(hg38_gwas_clean)[4] <- "N"

write.table(
    hg38_gwas_isotwas ,
    file = "PGC_UKB_23andMe_depression_genome-wide_isotwas_hg38_ver2.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
