## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library("data.table")
library("SummarizedExperiment")
library("here")
library("recount")
library("sva")
library("sessioninfo")
library("tidyverse")
library("getopt")

spec <- matrix(
    c('region', 'r', 1, 'character', 'Either Amygdala or SACC'),
    byrow = TRUE,
    ncol = 5
)
opt <- getopt(spec)

# opt$region <- "Amygdala"

if (tolower(opt$region) == "sacc") {
    opt$region = "sACC"
} else if (tolower(opt$region) == "amygdala" ||
           tolower(opt$region) == "amyg") {
    opt$region = "Amygdala"
} else {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Show the options used
message(paste(Sys.time(), "options used"))
print(opt)

#dir.create(here("hydeGoes_scSeq_mdd" , "MDD_vs_BP" , "isotwas" ),   paste0(opt$region, "_rda"), showWarnings = FALSE)
dir.create(paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , opt$region, "_rda"), showWarnings = FALSE)

## To avoid issues with running this code on qsub
data.table::setDTthreads(threads = 1)

## Find the samples for this project
load(here::here("hydeGoes_scSeq_mdd" ,  "MDD_vs_BP" , "data" , "rse_tx.Rdata"))

## 99 brains only have 1 sample either in one region or another
## 540 amyg samples
## 551 sACC samples
## 588 MDD
## 503 BP

# subset rse_tx to either Amygdala or sACC, whichever the user selected
rse_tx <- rse_tx[, colData(rse_tx)$BrainRegion == opt$region]

stopifnot(length(unique(rse_tx$BrNum)) == ncol(rse_tx))

# Load snpPCs

#/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/genotype_data/old
load(here::here("goesHyde_mdd_rnaseq" , "genotype_data", "old" ,  "goesHyde_bipolarMdd_Genotypes_mds.rda"))

# Read in original PLINK files
libd_bfile <-
    here::here("goesHyde_mdd_rnaseq" ,"genotype_data", "mdd_bpd", "maf01", "mdd_bpd_maf01.rsid")

## Read the LIBD fam data
libd_fam <- fread(
    paste0(libd_bfile, ".fam"),
    col.names = c(
        "famid",
        "w_famid",
        "w_famid_fa",
        "w_famid_mo",
        "sex_code",
        "phenotype"
    )
)

libd_fam$famKey <- paste0(libd_fam$famid, "_", libd_fam$w_famid)

## Filter the LIBD data to the one specific to this project
message(paste(Sys.time(), "processing", opt$region))
samp_file <- paste0("samples_to_extract_", opt$region, ".txt")

## Which samples have genotype data and MDS data?
samples_in_all <- intersect(rse_tx$genoSample, libd_fam$famKey)

################################
## Subset and save all key files
################################
rse_tx <-
    rse_tx[, rse_tx$genoSample %in% samples_in_all]

## Compute RPKM
#assays(rse_tx)$RPKM <- getRPKM(rse_tx, "Length")
#rse_tx_log = log2(assays(rse_tx)$tpm + 1)

## Compute Transcript PCs
message(Sys.time(), " computing Transcript PCs on log2(tpm + 1)")
pcaTranscript <- prcomp(t(log2(assays(rse_tx)$tpm + 1)))
save(pcaTranscript , file = paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , opt$region, "_rda/pcaTranscript.Rdata"))

message(Sys.time(), " determine how many Transcript PCs to adjust for")
mod <-
    model.matrix( ~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = colData(rse_tx))
kTranscript <- num.sv(log2(assays(rse_tx)$tpm + 1), mod)
stopifnot(kTranscript > 0)
TranscriptPCs <- pcaTranscript$x[, seq_len(kTranscript)]
save(TranscriptPCs, file = paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , opt$region, "_rda/TranscriptPCs.Rdata"))

## Add Transcript PCs to rse_tx
colData(rse_tx) <- cbind(colData(rse_tx), TranscriptPCs)

## Save for later
save(
    rse_tx,
    file = paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" ,
        opt$region,
        "_rda/",
        opt$region,
        "_hg38_rseTranscript_rawCounts_allSamples_n",
        ncol(rse_tx),
        ".Rdata"
    )
)

## Now extract the genotype data too
filter_m <- match(rse_tx$genoSample, libd_fam$famKey)
stopifnot(all(!is.na(filter_m)))
fwrite(libd_fam[filter_m, 1:2],
       ## can be more involved
       file =paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , samp_file),
       sep = "\t",
       col.names = FALSE)
newbfile_root <- paste0("mdd_bpd_maf01.rsid", opt$region)

dir.create(paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" , opt$region, "_duplicate_snps_bim"), showWarnings = FALSE)
newbfile <-
    here::here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , 
        paste0(opt$region, "_duplicate_snps_bim"),
        paste0(newbfile_root,
               "_duplicateSNPs")
    )

samp_file_address = paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/" ,  samp_file) 
## Extract
message(paste(Sys.time(), "running bfile extract for", newbfile))
system(
    paste(
        "plink --bfile",
        libd_bfile,
        "--keep",
        samp_file_address,
        "--make-bed --out",
        newbfile,
        " --memory 100000 --biallelic-only"
    )
)

## Check that we have the right data
newbfile_fam <- fread(
    paste0(newbfile, ".fam"),
    col.names = c(
        "famid",
        "w_famid",
        "w_famid_fa",
        "w_famid_mo",
        "sex_code",
        "phenotype"
    )
)

newbfile_fam$famKey <- paste0(newbfile_fam$famid, "_", newbfile_fam$w_famid)

check_m <-
    match(newbfile_fam$famKey,
          colData(rse_tx)$genoSample)

stopifnot(all(!is.na(check_m)))

## Re-run but now make the SNV names unique
dir.create(paste0("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas/",opt$region, "_unique_snps_bim"), showWarnings = FALSE)
newbfile_unique <-
      here::here("hydeGoes_scSeq_mdd/MDD_vs_BP/isotwas" , 
        paste0(opt$region, "_unique_snps_bim"),
        paste0(newbfile_root,
               "_uniqueSNPs")
    )

## Extract again (could also copy and rename, but it's fast this way)
message(paste(Sys.time(), "running bfile extract for", newbfile_unique))
system(
    paste(
        "plink --bfile",
        libd_bfile,
        "--keep",
        samp_file_address,
        "--make-bed --out",
        newbfile_unique,
        " --memory 100000 --biallelic-only"
    )
)


message(paste(Sys.time(), "reading the bim file", newbfile_unique))
bim <- fread(
    paste0(newbfile_unique, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

# > table(duplicated(bim$snp))
# 
# FALSE     TRUE 
# 10871656       10 

# > bim[duplicated(bim$snp),]$snp
# [1] "rs113644963"  "rs1209962010" "rs1344644938" "rs1165668794" "rs1236764134"
# [6] "rs1436504993" "rs1231536300" "rs542456559"  "rs74818883"   "rs755663741" 

## Make names unique
message(Sys.time(), " making the variant names unique")
bim$snp <- make.names(bim$snp, unique = TRUE)
stopifnot(all(!duplicated(bim$snp)))

## Ovewrite the PLINK bim file
fwrite(
    bim,
    file = paste0(newbfile_unique, ".bim"),
    sep = " ",
    col.names = FALSE
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

system("plink --version")
