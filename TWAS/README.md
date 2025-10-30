# What is FUSION-TWAS ?

- **Concept:** FUSION [[1](#references)] improved on original TWAS by better accounting for genetic architecture uncertainty and leveraging summary statistics exclusively.

- **Prediction Model:** FUSION addresses the uncertainty of genetic architecture (how many SNPs affect expression) by training **multiple distinct prediction models** (e.g., BLUP, BSLMM, LASSO, Elastic Net, Top-eQTL) for the same gene. It then selects the best-performing model (based on cross-validation $\text{R}^2$) for the TWAS test, or combines them (Omnibus test).

- **Association Step:** It is optimized for GWAS summary statistics (Z-scores/P-values) and is highly efficient. It uses these statistics combined with an LD reference panel and the selected expression weights ($\mathbf{w}$) to calculate the TWAS Z-score, effectively asking if the genetic component predicting expression is correlated with the trait.

- **Key Advantage:** Robustness across different genetic architectures due to testing multiple models.


# MDD vs. BP TWAS Workflow

This directory contains the computational workflow for the summary-TWAS (Transcriptome-Wide Association Study) conducted on MDD, BD, and control samples from the sACC (subgenual anterior cingulate cortex) and amygdala, as part of the broader MDD vs. BP (Bipolar Disorder) RNA-Seq project. 

Specifically, this workflow details the procedures used to generate TWAS results for both the amygdala and sACC tissues by leveraging the following public resources: the Bipolar Disorder GWAS data published by the Psychiatric Genomics Consortium (PGC) [[1](#references)], and the Bipolar Disorder GWAS data provided by 23andMe [[2](#references)].


## 1) Filter SNPs

```
sbatch 1_filter_snps_amygdala_gene.sh
# and/or
sbatch 1_filter_snps_sacc_gene.sh
```
The script first processes gene expression (RNA-seq) data, subsetting it to the chosen brain region, calculating RPKM (Reads Per Kilobase Million), and using Principal Component Analysis (PCA) and the sva package to estimate and correct for major confounding factors (Gene PCs), storing the processed expression data. Simultaneously, it uses the common samples to subset the raw genotype data (in PLINK format), extracts the relevant samples, and crucially, ensures that all Single Nucleotide Polymorphism (SNP) identifiers are unique by modifying and overwriting the PLINK .bim file. The final output is a set of carefully matched and corrected gene expression and genotype files (`.bed`, `.fam` and `.bim` files), ready for downstream TWAS modeling.

## 2) Build gene-level PLINK .bed files
```
sh 2_build_bims_Amygdala_gene.sh
# and/or
sh 2_build_bims_sACC_gene.sh
```
This script is the second critical phase of a Transcriptome-Wide Association Study (TWAS) pipeline and prepares the input files for the FUSION software. It first loads the previously processed gene expression data for a user-specified brain region, then applies a rigorous cleaning step using a linear model to regress out the effects of known and estimated confounders (like sex, population structure PCs, and expression PCs), resulting in **"cleaned expression" residuals**. Next, the script spatially filters the genes, keeping only those with at least **one SNP within a 1 Mb window**(500kb of coverage on each side of the gene). Finally, using **parallel processing (BiocParallel)**, it iterates through every remaining gene to create thousands of highly **localized PLINK file sets** (BED/BIM/FAM). In this crucial step, the gene's cleaned expression vector is saved as the phenotype in the .fam file, while the genotypes are subsetted to contain only the local SNPs, producing the required input format for building gene expression weight models.The outputs are in separate gene window subdirectories in `{subregion}_gene/bim_files/{subregion}_gene_{1:n gene windows}`.

## 3) Compute TWAS Weights individually

Here, we make use of Gusev et al's `FUSION.compute_weights.R` script [[2](#references)], which produces expression weights one gene at a time, taking in the above `build_bims.R` output. The flags are customizable and allow the user to specify multiple different models for testing. We set the heritability score P-value threshold to above `1` so heritability is not filtered. The output for this script can be found in `{subregion}_gene/out_files/gene_{1:n gene windows}`.

Because of limitation on SLURUM system, first by running bash script :
```
sh 3_0_preprocess_compute_weights_indv_Amygdala_full_gene_.sh 
and/or 
sh 3_0_preprocess_compute_weights_indv_sACC_full_gene_.sh
```
split the list of genes into smaller batches for parallel cluster submission, and create necessary symbolic links to define the input and output paths expected by the downstream FUSION TWAS software.

Then, using SLURM array scripts:
```
sbatch 3_compute_weights_indv_Amygdala_full_gene_1.sh to
sbatch 3_compute_weights_indv_Amygdala_full_gene_6.sh 
and/or 
sbatch 3_compute_weights_indv_sACC_full_gene_1.sh to
sbatch 3_compute_weights_indv_sACC_full_gene_6.sh
```

automates the parallel computation of gene expression weights for 4,500 individual genes in the Amygdala region by dynamically loading each gene's localized PLINK data and executing the FUSION TWAS R script with four different statistical models (Top1, BLUP, LASSO, Enet). 
## 4) Merge Individual TWAS gene weights
```
sbatch compute_weights_Amygdala_gene.sh
# and/or
sbatch compute_weights_sACC_gene.sh
```

The results from step #3 are merged together by generating an index `pos_info.Rdata`. `FUSION.profile_wgt.R` from Gusev et al.'s TWAS workflow is also called to produce a per-gene profile with summary of the data and model performance in `{subregion}_gene.profile.err`.

## 5) Process hg38 GWAS
```
sbatch run_process-gwas.sh
```

The GWAS file has already been converted to hg38, but requires some additional processing, such as the calculation of Z-score for each SNP and removal of unnecessary columns. Two histograms representing the frequency of Z-scores and log odds ratio are generated as well.

## 6) Apply Weights
```
sh apply_weights.sh
```

The weights generated in step #4 are used to perform expression imputation. `FUSION.assoc_test.R` is called, followed by `FUSION.post_process.R`, which conducts joint/conditional tests and produces corresponding plots. See more [here](http://gusevlab.org/projects/fusion/#typical-analysis-and-output) and [here](http://gusevlab.org/projects/fusion/#jointconditional-tests-and-plots).

## 7) TWAS results into a single Rdata file
```
sbatch run_read_twas_amygdala.sh
# and/or
sbatch run_read_twas_sacc.sh
```

All of the TWAS results from both the sACC and amygdala subregions are combined into a single Rdata file that is used for downstream analysis.

## 8) TWAS plot generation
```
sbatch run_generate_twas_plots.sh
```

This script generates a number of different plots and outputs their corresponding tables into `analysis/plots/` and `analysis/tables/`, respectively.

## References
1. Gusev, A., Ko, A., Shi, H., Bhatia, G., Chung, W., Penninx, B.W., Jansen, R., De Geus, E.J., Boomsma, D.I., Wright, F.A. and Sullivan, P.F., 2016. Integrative approaches for large-scale transcriptome-wide association studies. Nature genetics, 48(3), pp.245-252. https://www.nature.com/articles/ng.3506
2. O’Connell, K.S., Koromina, M., van der Veen, T., Boltz, T., David, F.S., Yang, J.M.K., Lin, K.H., Wang, X., Coleman, J.R., Mitchell, B.L. and McGrouther, C.C., 2025. Genomics yields biological and phenotypic insights into bipolar disorder. Nature, 639(8056), pp.968-975. https://www.nature.com/articles/s41586-024-08468-9  
3. Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118.
