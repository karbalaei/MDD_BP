# 🧬 What is isoTWAS?

- **The Concept**: Traditional TWAS evaluates associations at the aggregate gene level, which can obscure critical disease mechanisms. **isoTWAS** extends this framework by shifting the resolution to **isoform-specific expression (ISE)** and transcript-level features[[1](#references)]. This approach is uniquely powerful for distinguishing complex psychiatric conditions like MDD and BP, where alternative splicing and variable transcript usage often drive pathobiology even when total gene-level expression remains unchanged[[2](#references)].

- **Prediction Models** : Similar to gene-level frameworks, isoTWAS trains predictive models (such as Elastic Net, LASSO, and Blended architectures) to estimate the genetically regulated component of expression—but does so independently for individual transcripts or isoforms[[1](#references)].

- **Association Step**: By integrating transcript-specific predictive weights ($\mathbf{w}$) with major GWAS summary statistics (e.g., PGC or 23andMe) and an LD reference panel, isoTWAS computes isoform-level association $Z$-scores[[1](#references)]. This identifies specific transcript variations that are significantly correlated with phenotypic risk, allowing for high-resolution molecular mapping[[1 و 2](#references)].

## 🛠️ MDD vs. BP isoTWAS Workflow

Follow the sequential steps below to execute the transcript-level pipeline:.

## 1) Filter SNPs & Process Transcripts (`1_filter_snps_isotwast.R`)

Extracts and harmonizes the discovery transcriptomic and genotypic datasets for the chosen brain region (Amygdala or sACC). The script subsets the transcript-level expression data (`rse_tx`), filters for samples present in both the RNA-Seq and genotype cohorts, and calculates Transcript **Principal Components (PCs)** using log-transformed TPM data ($log_2(\text{tpm} + 1)$). The optimal number of hidden confounding factors to retain is determined dynamically via Surrogate Variable Analysis (`sva::num.sv`).Simultaneously, the script utilizes PLINK to isolate the corresponding sample genotypes under strict biallelic constraints (`--biallelic-only`). To protect downstream predictive modeling from pipeline crashes, it scans the genotypic data and systematically forces non-unique variant markers to become unique.

- **Outputs**: Cleaned transcript-level expression profiles (`rse_tx` containing estimated Transcript PCs) and matched, deduplicated PLINK binary genotype files (`.bed`, `.fam`, and `.bim`) formatted specifically for localized isoTWAS weight training.

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
1- Bhattacharya, A. et al. (2022) 'isoTWAS: Transcript-level transcriptome-wide association studies via structured lasso and blended learning', Genetics, 222(4), p.iyac145. doi:10.1093/genetics/iyac145.

2- Gandal, M. J. et al. (2018) 'Transcriptome-wide isoforms-level imbalances in the brains of individuals with severe mental illnesses', Science, 362(6420), p.eaat8127. doi:10.1126/science.aat8127.