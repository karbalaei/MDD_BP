library(sva)
# library("devtools")
# install_github("runehaubo/lmerTestR")
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(here)
library(purrr)
library(ggplot2)

rm(list = ls())


capabilities()

## load data

load(here::here( "data" , "rse_gene.Rdata"), verbose = TRUE)
load(here::here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)


## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Bipolar", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Bipolar", "MDD")]

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(rse_gene$PrimaryDx)

rse_gene$Dx<-relevel(rse_gene$Dx, "Bipolar")
rse_gene$Dx <- relevel(rse_gene$Dx, "Bipolar")

table(rse_gene$Dx)

# Bipolar     MDD
# 245     459


#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')

##### from run_wgcna_coombinedR script this was the model used for wgcna 
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr),
                        data=colData(rse_gene))

load(here::here("data","qSV_mat.Rdata"), verbose = TRUE)


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


save(modQsva, file = "rdas/modQsva_Michael_method.rda")

# [1] "(Intercept)"               "DxMDD"
#  [3] "BrainRegionsACC"           "AgeDeath"
#  [5] "DxMDD:BrainRegionsACC" "SexM"
#  [7] "snpPC1"                    "snpPC2"
#  [9] "snpPC3"                    "mitoRate"
# [11] "rRNA_rate"                 "totalAssignedGene"
# [13] "RIN"                       "ERCCsumLogErr"
# [15] "PC1"                       "PC2"
# [17] "PC3"                       "PC4"
# [19] "PC5"                       "PC6"
# [21] "PC7"                       "PC8"
# [23] "PC9"                       "PC10"
# [25] "PC11"                      "PC12"
# [27] "PC13"                      "PC14"
# [29] "PC15"                      "PC16"
# [31] "PC17"                      "PC18"
# [33] "PC19"                      "PC20"
# [35] "PC21"                      "PC22"


#########################
## load wgcna output ####
#########################

## load
load("rdas/constructed_network_signed_bicor_1.rda", verbose=TRUE)
# Loading objects:
#   net_list
#   net
#   fNames

#net is the main blockwide module command

### added code to create dendogram from WGCNA Tutorial
mergedColors = labels2colors(net$colors)

pdf(file = here::here("graphs" , "wgcna" ,  "dendrogram_Michael_method.pdf"))

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

dev.off() 
#########

# get colors LOOK AT WGCNA Instructions
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

dim(colorDat)
# [1] 18  4
### WILL CHANGE WITH AGE PROTECTED
colorDat

#18 4
# num          col Label numGenes
# ENSG00000227232.5    0         grey   ME0    14702
# ENSG00000228794.8    1    turquoise   ME1     2710
# ENSG00000225630.1    2         blue   ME2     2501
# ENSG00000162572.19   3        brown   ME3     2333
# ENSG00000142583.17   4       yellow   ME4      680
# ENSG00000176022.4    5        green   ME5      409
# ENSG00000142609.17   6          red   ME6      327
# ENSG00000048707.13   7        black   ME7      278
# ENSG00000187730.7    8         pink   ME8      211
# ENSG00000069424.14   9      magenta   ME9      155
# ENSG00000107404.18  10       purple  ME10      145
# ENSG00000162545.5   11  greenyellow  ME11      144
# ENSG00000162576.16  12          tan  ME12      137
# ENSG00000077254.14  13       salmon  ME13      132
# ENSG00000067606.16  14         cyan  ME14      127
# ENSG00000162585.16  15 midnightblue  ME15      111
# ENSG00000142627.12  16    lightcyan  ME16       62
# ENSG00000183888.4   17       grey60  ME17       48



######################
## gene ontology
### split genes by each color 
gList = split(rowData(rse_gene)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene)$EntrezID
univ = as.character(univ[!is.na(univ)])


ont <- c("BP", "CC", "MF" , "ALL")
names(ont) <- ont


plan <- tidyr::expand_grid(
  gene_list = gList,
  ontology = ont
)


Enrichment_results <-  pmap(plan ,  ~enrichGO(.x,
                                              universe = univ, OrgDb = "org.Hs.eg.db",
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


walk2(names(all_plots), all_plots, ~ggsave(filename = here::here("graphs" , "wgcna" , "enrichment" ,  paste0(.x , "_Michael_method.jpg")), plot = .y, 
                                           height = 14, width = 9))



pdf(here::here("graphs" , "wgcna" , "enrichment" , "Enrichment results  WGNCA MDD vs. BP_Michael_method.pdf") , height = 14 , width = 9)

map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results,file = here::here("results", "wgcna_enrichGO_Michael_method.rda"))



##############
## associate eigengenes with brain region
m = modQsva[,1:5] # this is what was protected
colnames(m)
#### look up MEs in WGCNA 
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col] #change order of columns
dim(MEs)
## check## lmer intercept -- random effect -- check what -1 means 
statList = lapply(MEs, function(x) summary(lmer(x ~ m + (1|rse_gene$BrNum) - 1))$coef)

statList[[1]]


# Estimate   Std. Error       df    t value   Pr(>|t|)
# m(Intercept)           -0.0011858385 0.0060608647 438.8080 -0.1956550 0.84497077
# mDxMDD                 -0.0045760217 0.0041498116 592.4354 -1.1027059 0.27060265
# mBrainRegionsACC        0.0051267527 0.0034255550 360.8223  1.4966196 0.13536646
# mAgeDeath               0.0001061846 0.0001178341 374.1109  0.9011367 0.36809541
# mDxMDD:BrainRegionsACC -0.0087350405 0.0042169799 354.4810 -2.0713972 0.03904471



# modified extract information from lme object
MDDEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
regionEffect = as.data.frame(t(sapply(statList, function(x) x[3,])))
## now protect age we cannot 
ageEffect = as.data.frame(t(sapply(statList, function(x) x[4,])))  
#interaction effect
intEffect = as.data.frame(t(sapply(statList, function(x) x[5,])))
#colnames(MDDEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
#	"slope", "se", "df", "t", "pvalue")

colnames(MDDEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
	"slope", "se", "df", "t", "pvalue")

print_effect <- function(x) { signif(cbind(x, FDR = p.adjust(x$pvalue, 'fdr')), 3) }
print_effect(MDDEffect)


                # slope      se  df      t   pvalue      FDR
# grey         -0.00458 0.00415 592 -1.100 2.71e-01 3.75e-01
# turquoise     0.00165 0.00191 682  0.866 3.87e-01 4.64e-01
# blue          0.00532 0.00218 693  2.440 1.51e-02 3.40e-02
# brown        -0.02370 0.00358 682 -6.630 7.01e-11 6.31e-10
# yellow        0.00767 0.00414 602  1.850 6.42e-02 1.05e-01
# green         0.00538 0.00409 691  1.310 1.89e-01 2.84e-01
# red           0.01810 0.00395 688  4.580 5.60e-06 3.36e-05
# black        -0.01440 0.00407 695 -3.530 4.46e-04 1.61e-03
# pink          0.02670 0.00393 695  6.780 2.51e-11 4.51e-10
# magenta       0.00141 0.00415 666  0.339 7.34e-01 7.78e-01
# purple       -0.00877 0.00408 663 -2.150 3.19e-02 6.38e-02
# greenyellow   0.00176 0.00422 696  0.417 6.77e-01 7.62e-01
# tan          -0.00109 0.00417 683 -0.262 7.93e-01 7.93e-01
# salmon       -0.01320 0.00359 697 -3.690 2.46e-04 1.11e-03
# cyan          0.01280 0.00396 693  3.220 1.32e-03 3.43e-03
# midnightblue -0.01260 0.00389 536 -3.230 1.33e-03 3.43e-03
# lightcyan    -0.00832 0.00396 499 -2.100 3.59e-02 6.46e-02
# grey60       -0.00352 0.00398 600 -0.885 3.76e-01 4.64e-01




#signif(ageEffect, 3)
print_effect(ageEffect)
              # slope      se  df      t   pvalue      FDR
# grey          1.06e-04 1.18e-04 374  0.9010 0.368000 0.42500
# turquoise     1.68e-04 4.81e-05 362  3.5000 0.000523 0.00235
# blue         -4.70e-05 5.33e-05 361 -0.8820 0.378000 0.42500
# brown         3.30e-04 9.00e-05 361  3.6600 0.000285 0.00171
# yellow        3.48e-04 1.16e-04 373  2.9900 0.002970 0.00892
# green        -5.34e-06 1.00e-04 353 -0.0532 0.958000 0.95800
# red           3.79e-04 9.85e-05 376  3.8500 0.000137 0.00171
# black         2.31e-04 9.86e-05 365  2.3400 0.019700 0.04280
# pink         -1.44e-04 9.48e-05 357 -1.5200 0.130000 0.18000
# magenta      -2.48e-04 1.07e-04 338 -2.3100 0.021400 0.04280
# purple       -1.02e-04 1.06e-04 350 -0.9580 0.338000 0.42500
# greenyellow   1.16e-05 1.01e-04 342  0.1150 0.909000 0.95800
# tan           3.93e-04 1.05e-04 370  3.7400 0.000209 0.00171
# salmon        1.58e-04 8.57e-05 365  1.8500 0.065500 0.10700
# cyan         -1.48e-04 9.59e-05 323 -1.5400 0.125000 0.18000
# midnightblue -2.38e-04 1.16e-04 371 -2.0600 0.040200 0.07240
# lightcyan    -4.09e-04 1.21e-04 372 -3.3900 0.000769 0.00277
# grey60       -2.72e-04 1.12e-04 376 -2.4200 0.016000 0.04100


print_effect(intEffect)
#           # slope      se  df      t   pvalue      FDR
# grey         -0.008740 0.00422 354 -2.0700 0.0390 0.316
# turquoise    -0.006030 0.00244 362 -2.4700 0.0140 0.251
# blue         -0.000172 0.00292 365 -0.0588 0.9530 0.953
# brown        -0.005190 0.00459 361 -1.1300 0.2590 0.777
# yellow       -0.003990 0.00430 355 -0.9270 0.3540 0.795
# green        -0.003200 0.00543 357 -0.5890 0.5560 0.795
# red          -0.010000 0.00515 378 -1.9400 0.0526 0.316
# black        -0.002400 0.00550 370 -0.4360 0.6630 0.795
# pink         -0.005220 0.00533 363 -0.9790 0.3280 0.795
# magenta       0.003600 0.00513 334  0.7010 0.4840 0.795
# purple       -0.002660 0.00499 345 -0.5330 0.5940 0.795
# greenyellow   0.002610 0.00577 349  0.4530 0.6510 0.795
# tan          -0.007440 0.00535 371 -1.3900 0.1650 0.741
# salmon       -0.002900 0.00492 372 -0.5900 0.5560 0.795
# cyan         -0.006270 0.00534 329 -1.1700 0.2410 0.777
# midnightblue -0.001210 0.00342 343 -0.3530 0.7240 0.815



##################
# make boxplots ##
lab = paste0(substr(rse_gene$BrainRegion,1,4), ":", rse_gene$Dx)
table(lab)
# Amyg:Bipolar     Amyg:MDD sACC:Bipolar     sACC:MDD
#      122          231          123          228

### changed BP to MDD labels
lab = factor(lab, levels = c("sACC:Bipolar", "sACC:MDD", "Amyg:Bipolar", "Amyg:MDD"))

pdf(here::here("graphs" , "wgcna" , "MEs_vs_dx_Michael_method.pdf"),w=8,h=6, useDingbats = FALSE)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=0.8,cex.lab=1,cex.main = 1.4)
for(i in 1:ncol(MEs)) {
	boxplot(MEs[,i] ~ lab, outline = FALSE, xlab="",
		ylim = quantile(unlist(MEs),c(0.001,0.999)),main = colnames(MEs)[i],
		names = gsub(":", "\n", levels(lab)), ylab = "Module Eigengene")
	points(MEs[,i] ~ jitter(as.numeric(lab),amount=0.1), pch=21, bg=lab)
	legend("top", c(paste0("Region p=", signif(regionEffect[i,5],3)), 
		paste0("Dx p=", signif(MDDEffect[i,5],3))),cex=1)
}
dev.off()

colnames(m)
# [1] "(Intercept)"               "DxControl"                
# [3] "BrainRegionsACC"           "AgeDeath"                 
# [5] "DxControl:BrainRegionsACC"

clean_ME = t(cleaningY(t(MEs), m, P=2))
pdf(here::here("graphs" , "wgcna" ,"clean_MEs_vs_dx_Michael_method.pdf"),w=5,h=6, useDingbats = FALSE)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=0.8,cex.lab=1,cex.main = 1.4)
for(i in 1:ncol(clean_ME)) {
	boxplot(clean_ME[,i] ~ rse_gene$Dx, outline = FALSE, xlab="",
		ylim = quantile(unlist(clean_ME),c(0.001,0.999)),main = colnames(MEs)[i],
		ylab = "Module Eigengene (Adj)")
	points(clean_ME[,i] ~ jitter(as.numeric(rse_gene$Dx),amount=0.1), pch=21, bg=lab)
	legend("top", paste0("Dx p=", signif(MDDEffect[i,5],3)),cex=1)
}
dev.off()
	


### enrichment of DEGs #####

load(here::here("Final results" , "Raw R file results" ,  "BP_MDD_All_results_Micheal_method_Final.RDS"))

# load("../case_control/interaction_model_results.rda")
# ## remove everything but gene level
# rm(outExon_bothRegion, outJxn_bothRegion, outTx_bothRegion)

### Gene ####
#### All #####

#Amygdala

identical(rownames(Outlist$gene$amyg), fNames) # TRUE
tt = table(net$colorsLab, Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05)

tt = tt[rownames(MDDEffect),]

pdf(here::here ("Final results" , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect Amygdala_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
	tab = table(net$colorsLab == cc,Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05)
	c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect Amygdala_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_Amygdala_table_Michael_method.csv" ) , sep = "," )

#sACC

identical(rownames(Outlist$gene$sacc), fNames) # TRUE
tt = table(net$colorsLab, Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05)
tt = tt[rownames(MDDEffect),]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect sACC_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net$colorsLab == cc,Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05)
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect sACC_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_sACC_table_Michael_method.csv" ) , sep = "," )



##### Up #####

#Amygdala

identical(rownames(Outlist$gene$amyg), fNames) # TRUE
tt = table(net$colorsLab, c(Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05 & Outlist$gene$amyg$PrimaryDxBipolar > 0 ))
tt = tt[rownames(MDDEffect),]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect Amygdala_UP_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net$colorsLab == cc, c(Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05 & Outlist$gene$amyg$PrimaryDxBipolar > 0 ))
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect Amygdala_UP_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_Amygdala_table_UP_Michael_method.csv" ) , sep = "," )

#sACC

identical(rownames(Outlist$gene$sacc), fNames) # TRUE
tt = table(net$colorsLab, c(Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05 & Outlist$gene$sacc$PrimaryDxBipolar > 0 ))
tt = tt[rownames(MDDEffect),]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect sACC_UP_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net$colorsLab == cc, c(Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05 & Outlist$gene$sacc$PrimaryDxBipolar > 0 ))
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect sACC_UP_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_sACC_table_UP_Michael_method.csv" ) , sep = "," )


##### DOWN ######

#Amygdala

identical(rownames(Outlist$gene$amyg), fNames) # TRUE
tt = table(net$colorsLab, c(Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05 & Outlist$gene$amyg$PrimaryDxBipolar < 0 ))
tt = tt[rownames(MDDEffect),]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect Amygdala_DOWN_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net$colorsLab == cc, c(Outlist$gene$amyg$q_PrimaryDxBipolar < 0.05 & Outlist$gene$amyg$PrimaryDxBipolar < 0 ))
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect Amygdala_DOWN_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_Amygdala_table_DOWN_Michael_method.csv" ) , sep = "," )

#sACC

identical(rownames(Outlist$gene$sacc), fNames) # TRUE
tt = table(net$colorsLab, c(Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05 & Outlist$gene$sacc$PrimaryDxBipolar < 0 ))
tt = tt[rownames(MDDEffect),]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "prop.table vs MDDEffect sACC_DOWN_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

dev.off()
# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net$colorsLab == cc, c(Outlist$gene$sacc$q_PrimaryDxBipolar < 0.05 & Outlist$gene$sacc$PrimaryDxBipolar < 0 ))
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

pdf(here::here("Final results"  , "WGCNA" , "MDD effect graphs" , "-log10Pvalue chi-sq vs MDDEffect sACC_DOWN_Michael_method.pdf"),w=9,h=6, useDingbats = FALSE)

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

dev.off()

write.table(deModEnrich , file = here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,  "deModEnrich_sACC_table_DOWN_Michael_method.csv" ) , sep = "," )






##### oTHER FEtures ####

identical(rownames(Outlist$gene$amyg), fNames)

length(rownames(Outlist$gene$amyg))  #25212
length(fNames) #25212

sum(fNames %in% rownames(Outlist$gene$amyg)  ) #24784


sum(rownames(Outlist$gene$amyg) %in% fNames) #24784

net$fNames <- fNames


#### Genes #####

#name = "Gene_Amygdala" # i =1 j = 1
 name = "Gene_sACC" # i =1 j = 2
# 
i =1
# 
j=1

rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)

Enrichment_fun_gene <- function( i , j , name) {
  
  Avi_id <-   which(net$fNames %in% rownames(Outlist[[i]][[j]]))
  
  net_new <-  list()
  
  net_new$colorsLab <-  net$colorsLab[Avi_id]
  
  tt = table(net_new$colorsLab, Outlist[[i]][[j]]$q_PrimaryDxBipolar < 0.05)
  tt = tt[rownames(MDDEffect),]
  
  prop.table_tt  = prop.table(tt,1)
  
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
  
  
  # manual chi-sq?
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_new$colorsLab == cc,Outlist[[i]][[j]]$q_PrimaryDxBipolar < 0.05)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  numGenes <-  table(net_new$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(where(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  
  write.csv(prop.table_tt , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs",  paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  
  plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs",  paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs",  paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
  
}

gene_plan <-  data.frame( i = 1 , j = c(2 , 1) , name = c("Gene_Amygdala" , "Gene_sACC"))

pmap(gene_plan , Enrichment_fun_gene)

#### Exon #####

rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)

name = "Exon_Amygdala" # i =2 j = 1
name = "Exon_sACC" # i =1 j = 2

i =2

j=1

Enrichment_fun_exon <- function( i , j , name) {
  
  common_id <-   intersect(net$fNames , Outlist[[i]][[j]]$common_gene_id)
  
  Avi_id <-   which(net$fNames %in%  common_id)
  
  Outlist_common <-  Outlist[[i]][[j]][Outlist[[i]][[j]]$common_gene_id %in% common_id ,]
  
  Outlist_DEf <-  Outlist_common[Outlist_common$q_PrimaryDxBipolar < 0.05 , 'common_gene_id' ] %>% unique()
  
  univers <-  Outlist_common$common_gene_id %>% unique() %>% as.data.frame() %>%
    setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Outlist_DEf , TRUE , FALSE))
  
  net_new <-  list()
  
  net_new$colorsLab <-  net$colorsLab[Avi_id]
  
  tt = table(net_new$colorsLab, univers$DEf )
  tt = tt[rownames(MDDEffect),]
  
  prop.table_tt  = prop.table(tt,1)
  
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
  
  
  # manual chi-sq?
  
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_new$colorsLab == cc,univers$DEf)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  numGenes <-  table(net_new$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(here(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  
  
  write.csv(prop.table_tt , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
}

exon_plan <-  data.frame( i = 2 , j = c(2 , 1) , name = c("Exon_Amygdala" , "Exon_sACC"))

pmap(exon_plan , Enrichment_fun_exon)

#### Junction #####

rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)

name = "Junction_Amygdala" # i =2 j = 1
name = "Junction_sACC" # i =1 j = 2

i =3

j=1

Enrichment_fun_junction <- function( i , j , name) {
  
  common_id <-   intersect(net$fNames , unique(Outlist[[i]][[j]]$common_gene_id))
  
  Avi_id <-   which(net$fNames %in%  common_id)
  
  Outlist_common <-  Outlist[[i]][[j]][Outlist[[i]][[j]]$common_gene_id %in% common_id ,]
  
  Outlist_DEf <-  Outlist_common[Outlist_common$q_PrimaryDxBipolar < 0.05 , 'common_gene_id' ] %>% unique()
  
  
  univers <-  Outlist_common$common_gene_id %>% unique() %>% as.data.frame() %>%
    setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Outlist_DEf , TRUE , FALSE)) 
  net_new <-  list()
  
  net_new$colorsLab <-  net$colorsLab[Avi_id]
  
  
  
  #Avi_id2 <-   which(net$fNames %in% univers)
  
  tt = table(net_new$colorsLab, univers$DEf )
  tt = tt[rownames(MDDEffect),]
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  prop.table_tt  = prop.table(tt,1)
  
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num)  %>% mutate(across(where(is.numeric), round, 3)) 
  
  # manual chi-sq?
  
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_new$colorsLab == cc,univers$DEf)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  numGenes <-  table(net_new$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(where(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  write.csv(prop.table_tt , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
}

junction_plan <-  data.frame( i = 3 , j = c(2 ,1) , name = c("Junction_Amygdala" , "Junction_sACC"))

pmap(junction_plan , Enrichment_fun_junction)

#### Transcript #####

# name = "Junction_Amygdala" # i =2 j = 1
# name = "Junction_sACC" # i =1 j = 2
# 
# i =4
# 
# j=1

rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)


Enrichment_fun_transcript <- function( i , j , name) {
  
  common_id <-   intersect(net$fNames , unique(Outlist[[i]][[j]]$common_gene_id))
  
  Avi_id <-   which(net$fNames %in%  common_id)
  
  Outlist_common <-  Outlist[[i]][[j]][Outlist[[i]][[j]]$common_gene_id %in% common_id ,]
  
  Outlist_DEf <-  Outlist_common[Outlist_common$q_PrimaryDxBipolar < 0.05 , 'common_gene_id' ] %>% unique()
  
  
  univers <-  Outlist_common$common_gene_id %>% unique() %>% as.data.frame() %>%
    setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Outlist_DEf , TRUE , FALSE))
  net_new <-  list()
  
  net_new$colorsLab <-  net$colorsLab[Avi_id]
  
  
  
  #Avi_id2 <-   which(net$fNames %in% univers)
  
  tt = table(net_new$colorsLab, univers$DEf )
  tt = tt[rownames(MDDEffect),]
  
  prop.table_tt  = prop.table(tt,1)
  
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
  
  
  # manual chi-sq?
  
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_new$colorsLab == cc,univers$DEf)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  
  numGenes <-  table(net_new$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(where(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  
  write.csv(prop.table_tt , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("Final results"  , "WGCNA" , "MDD effect graphs" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("Final results"  , "WGCNA" , "MDD effect graphs" ,    paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
}

transcript_plan <-  data.frame( i = 4 , j = c(1 , 2) , name = c("Transcript_Amygdala" , "Transcript_sACC"))

pmap(transcript_plan , Enrichment_fun_transcript)



# 
# ## Reproducibility information
 print('Reproducibility information:')
 Sys.time()
 proc.time()
 options(width = 120)
 session_info()

