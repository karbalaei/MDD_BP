library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(gwascat)

hg19_gwas <-
    fread(
        "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/hg19_BP/bip2024_multianc_no23andMe.txt"
    )

glimpse(hg19_gwas)


# Add 'chr' prefix 
hg19_gwas$CHR <- ifelse(!grepl("^chr", hg19_gwas$CHR), paste0("chr", hg19_gwas$CHR), hg19_gwas$CHR)

glimpse(hg19_gwas)

#Convert your data to a GRanges object
granges_hg19 <- GRanges(
    seqnames = hg19_gwas$CHR,
    ranges = IRanges(start = hg19_gwas$BP, end = hg19_gwas$BP),
    snp_id = hg19_gwas$SNP
)

#Download the chain file and import it

setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/hg19_BP")
chain_file <- "hg19ToHg38.over.chain"

if (file.exists(chain_file)) {

  chain <- import.chain(chain_file)

} else {

download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", destfile = paste0(chain_file, ".gz"))
R.utils::gunzip(paste0(chain_file, ".gz"))

chain <- import.chain(chain_file)

}



#Perform the liftOver

seqlevelsStyle(granges_hg19) = "UCSC"  # necessary


granges_hg38 <- liftOver(granges_hg19, chain)

granges_hg38 <- unlist(granges_hg38)

#convert hg38 to dataframe

hg38 <- data.frame(
    CHR = seqnames(granges_hg38),
    BP = start(granges_hg38),
    SNP = granges_hg38$snp_id
)


hg38_gwas = hg19_gwas[ , c( 1, 4:9)]  %>% merge(hg38 , by = "SNP") 


# read bim file with unique rsIDs


uniq_bim <-
    fread(
        here::here("twas", "Amygdala_unique_snps_bim", "mdd_bpd_maf01.rsidAmygdala_uniqueSNPs.bim")
    )


colnames(uniq_bim) <- c("CHR", "SNP", "dummy", "BP", "A1", "A2")

# check for overlap between gwas and unique bim file
uniq_bim$new_test <- paste0(uniq_bim$CHR, "_", uniq_bim$BP)

hg38_gwas$new_test <-
    paste0(gsub("chr", "" , hg38_gwas$CHR), "_", hg38_gwas$BP)

table(hg38_gwas$new_test %in% uniq_bim$new_test)
#  FALSE    TRUE
# 46938 6310719


hg38_gwas <-
    merge(hg38_gwas[,-"SNP"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas$new_test <- NULL

col_order <- c(
    "CHR",
    "SNP",
    "BP",
    "A1",
    "A2",
    "INFO",
    "OR",
    "SE",
    "P")

hg38_gwas <- hg38_gwas[, ..col_order]



# calculate z score

hg38_gwas$Z <- log(hg38_gwas$OR) / hg38_gwas$SE



setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/MDD_vs_BP/twas/hg38")
pdf(file = "bip2024_multianc_no23andMe_HIST_hg38.pdf", useDingbats = FALSE)

hist(log(hg38_gwas$OR), col = "gold" , main = "Histogram of natural logarithm of OR", xlab = "Natural logarithm of OR")
hist(hg38_gwas$Z, col = "darkorange" , main = "Histogram of Z", xlab = "Z")

dev.off()

hg38_gwas_clean <-
    hg38_gwas[, c("SNP", "A1", "A2", "P", "Z")]
# colnames(hg38_gwas_clean)[4] <- "N"

write.table(
    hg38_gwas_clean,
    file = "bip2024_multianc_no23andMe_CLEAN_hg38.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

