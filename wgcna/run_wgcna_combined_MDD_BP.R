#######
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)
library(sessioninfo)
library(here)
library(parallel)
library(doParallel)

#qrsh -l mem_free=100G,h_vmem=100G

## multithread
#setallowWGCNAThreads(8)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(40)


#dir.create("rdas", showWarnings = FALSE)

## load data
#load(here("exprs_cutoff","rse_gene.Rdata"))

load(here( "data" , "rse_gene.Rdata"), verbose = TRUE)
load(here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)
load(here("data","qSV_mat.Rdata"), verbose = TRUE)

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 459     387     245

table(cov_rse$PrimaryDx)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Bipolar", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Bipolar", "MDD")]
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 463     387       0

qSV_mat


rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(cov_rse$PrimaryDx)

rse_gene$Dx<-relevel(rse_gene$Dx, "Bipolar")
rse_gene$Dx <- relevel(rse_gene$Dx, "Bipolar")

table(rse_gene$Dx)

#Bipolar     MDD 
#    245     459 


## add ancestry
#no need -- already in new rse_object
#load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)
#class(mds)
#dim(mds)
#corner(mds)
#table(rse_gene$BrNum %in% rownames(mds))

# FALSE  TRUE
# 12   836
## 6 individual, 12 brain regions absent in mds file

#table(rse_gene$BrNum %in% rownames(mds), rse_gene$BrainRegion)
# Amygdala sACC
# FALSE        6    6
# TRUE       415  421

## keep samples with genotypes
#rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
#cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

dim(rse_gene)
# [1] 25212   704
addmargins(table(rse_gene$Dx, rse_gene$BrainRegion))

#         Amygdala sACC Sum
# Bipolar      122  123 245
# MDD          231  228 459
# Sum          353  351 704


#### reorders according to rse_gene$BrNum
#mds = mds[rse_gene$BrNum,1:5]
#dim(mds)
# [1] 836   5

#colnames(mds) = paste0("snpPC", 1:5)
#colData(rse_gene) = cbind(colData(rse_gene), mds)
#colData(cov_rse) = cbind(colData(cov_rse), mds)

###########
#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')


##########
## model #
##########

### came from qSV_model_DE_analysis.R
# modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate +
                            # totalAssignedGene + RIN, data = colData(rse_gene))


### removed ERCCsumLogErr term
### REPLACED ERCCsumLogErr term

# modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
# 	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
# 	data=colData(rse_gene))

#update to include abs ERCC

modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr),
	data=colData(rse_gene))

load(here("data","qSV_mat.Rdata"), verbose = TRUE)


### counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
print(k) # 22
qSV_mat = prcomp(t(degExprs))$x[,1:k]

## join and move around region, dx and interaction for cleaning

colnames(modJoint)
# [1] "(Intercept)"           "DxMDD"                 "BrainRegionsACC"
# [4] "AgeDeath"              "SexM"                  "snpPC1"
# [7] "snpPC2"                "snpPC3"                "snpPC4"
# [10] "snpPC5"                "mitoRate"              "rRNA_rate"
# [13] "totalAssignedGene"     "RIN"                   "abs(ERCCsumLogErr)"
# [16] "DxMDD:BrainRegionsACC"


modQsva = cbind(modJoint[,c(1:4,16,5:15)], qSV_mat)

colnames(modQsva)


# [1] "(Intercept)"           "DxMDD"                 "BrainRegionsACC"
# [4] "AgeDeath"              "DxMDD:BrainRegionsACC" "SexM"
# [7] "snpPC1"                "snpPC2"                "snpPC3"
# [10] "snpPC4"                "snpPC5"                "mitoRate"
# [13] "rRNA_rate"             "totalAssignedGene"     "RIN"
# [16] "abs(ERCCsumLogErr)"    "PC1"                   "PC2"
# [19] "PC3"                   "PC4"                   "PC5"
# [22] "PC6"                   "PC7"                   "PC8"
# [25] "PC9"                   "PC10"                  "PC11"
# [28] "PC12"                  "PC13"                  "PC14"
# [31] "PC15"                  "PC16"                  "PC17"
# [34] "PC18"                  "PC19"                  "PC20"
# [37] "PC21"                  "PC22"

## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)

### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)

#geneExprsClean_sex = cleaningY(geneExprs, modQsva, P=6)

#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,
                               networkType = "signed", verbose = 5)
sftthresh1$powerEstimate
#10
save(sftthresh1, file =  here("rdas/power_object.rda"))

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase =here("rdas/wgcna_signed_TOM"))
fNames = rownames(geneExprs)

### saved in later command
save(net, fNames, file = here("rdas/constructed_network_signed_bicor_1.rda"))

########################
## by region remove 
modRegion =  model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr),
	data=colData(rse_gene))

colnames(modRegion)
# [1] "(Intercept)"       "DxControl"         "AgeDeath"
# [4] "SexM"              "snpPC1"            "snpPC2"
# [7] "snpPC3"            "mitoRate"          "rRNA_rate"
# [10] "totalAssignedGene" "RIN"

###jaffe lab function  splits by brain region 
rIndexes = splitit(rse_gene$BrainRegion)

### output of splitit 2 list 
## clean by region ## changed P=3 to P=2 since not clear why to protect AgeDeath
### loop for each brain region


geneExprs_list = mclapply(rIndexes, function(ii) {
##extract the deg exp for each brain region samples
  	degExprs = log2(assays(cov_rse[,ii])$count+1)
#run number of svs for each brain region samples  	
	k = num.sv(degExprs, modRegion[ii,])
#combine model terms/matrix with qSVs *(prcom is a principal component function)l PC on degradation data 	
	m = cbind(modRegion[ii,], prcomp(t(degExprs))$x[,1:k])
#regress out all the effects for gene expression except for intercept and DxControl for each brain region  	
	cleaningY(geneExprs[,ii], m, P=2)
},mc.cores=2)

## threshold
thresh_list = mclapply(geneExprs_list, function(y) {
	pickSoftThreshold(t(y), powerVector = powers,
                               networkType = "signed", verbose = 5)
},mc.cores=2)

## networks
net_list = lapply(1:2, function(i) {
  blockwiseModules(t(geneExprs_list[[i]]), power = thresh_list[[i]]$powerEstimate,
                   networkType = "signed", minModuleSize = 30,corType="bicor",
                   reassignThreshold = 0, mergeCutHeight = 0.25,
                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                   saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                   saveTOMFileBase = here("rdas" , paste0("wgcna_signed_TOM_region",
								names(rIndexes)[i])))
})

save(net_list, net, fNames, file = here("rdas/constructed_network_signed_bicor_2.rda"))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
