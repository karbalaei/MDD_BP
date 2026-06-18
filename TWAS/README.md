# Transcriptome-Wide Association Studies (TWAS) Pipeline

This directory contains the computational workflow for the summary-based Transcriptome-Wide Association Study (TWAS) conducted on Bipolar Disorder (BP) and Major Depressive Disorder (MDD) cohorts. This analysis utilizes RNA-Seq expression profiles from two critical brain regions—the **amygdala** and the **subgenual anterior cingulate cortex (sgACC)**—to identify genetically regulated gene expression variations associated with these mood disorders.

---

## 🧬 What is FUSION-TWAS?

* **The Concept:** FUSION [[1](#references)] improves upon traditional TWAS frameworks by better accounting for the uncertainty of complex local genetic architectures and leveraging summary statistics exclusively.
* **Prediction Models:** To capture varying underlying genetic models (e.g., sparse vs. highly polygenic effects), FUSION trains **multiple distinct prediction models** simultaneously for each gene:
  * **BLUP** (Best Linear Unbiased Predictor)
  * **BSLMM** (Bayesian Sparse Linear Mixed Model)
  * **LASSO** / **Elastic Net**
  * **Top-eQTL** (single top SNP model)
  
  The framework uses cross-validation $R^2$ to select the best-performing model for the association test, or combines them via an Omnibus test.
* **Association Step:** FUSION is highly optimized for GWAS summary statistics ($Z$-scores or $P$-values). It combines these stats with a linkage disequilibrium (LD) reference panel and pre-computed expression weights ($\mathbf{w}$) to calculate a final TWAS $Z$-score, assessing whether the predicted genetic component of gene expression correlates with the trait.

---

## 🛠️ MDD vs. BP TWAS Workflow

This workflow processes local RNA-Seq datasets alongside public-facing genetic datasets, leveraging the Bipolar Disorder and MDD GWAS summary data from the Psychiatric Genomics Consortium (PGC) and 23andMe meta-analyses [[2](#references)].

Follow the sequential steps below to execute the pipeline:

### 1) Filter SNPs & Match Formats
Processes the raw genotype and discovery RNA-Seq datasets. It subsets data to the specific brain region, computes RPKM, and corrects for latent confounding factors using Principal Component Analysis (PCA) and Surrogate Variable Analysis via the `sva` package. Concurrently, it extracts corresponding samples from the raw PLINK files and standardizes SNP identifiers in the `.bim` file.
* **Outputs:** Matched, cleaned `.bed`, `.fam`, and `.bim` files ready for expression modeling.

```bash
# Execute for the desired brain region:
sbatch 1_filter_snps_amygdala_gene.sh
# and/or
sbatch 1_filter_snps_sacc_gene.sh
```


## 2) Build Gene-Level PLINK Files
Prepares localized window boundaries for FUSION. This script uses a linear model to regress out known covariates (sex, population structure PCs, expression PCs) to yield "cleaned expression" residuals. It then applies a structural spatial filter, capturing all SNPs within a 1 Mb cis-regulatory window (500 kb upstream/downstream of the gene body) using parallel processing (BiocParallel).

- **Outputs** : Thousands of independent gene subdirectories stored in {subregion}_gene/bim_files/ containing localized PLINK inputs with expression residuals designated as the phenotype.

```
sh 2_build_bims_Amygdala_gene.sh
# and/or
sh 2_build_bims_sACC_gene.sh
```

## 3) Compute Local Expression Weights

Executes the core FUSION weight computation pipeline (FUSION.compute_weights.R) to model cis-eQTL genetic architectures [1]. The heritability $P$-value threshold is set to 1 to prevent aggressive filtering of low-heritability features at early stages.

- **Note on Cluster Scheduling** : Due to SLURM step limitations, individual batch arrays are dynamically chunked to execute smoothly.
```
# Step 3.0: Preprocess and split the gene catalog into batches
sh 3_0_preprocess_compute_weights_indv_Amygdala_full_gene_.sh 
sh 3_0_preprocess_compute_weights_indv_sACC_full_gene_.sh

# Step 3.1: Submit SLURM job arrays across the generated blocks
sbatch 3_compute_weights_indv_Amygdala_full_gene_1.sh # through 6
sbatch 3_compute_weights_indv_sACC_full_gene_1.sh     # through 6

```

## 4) Merge Individual Gene Weights
Aggregates the individual expression weights calculated across the thousands of parallel cluster slots. It references the underlying architecture by constructing a global index file (pos_info.Rdata) and generates model fit diagnostics via FUSION performance profiles.

- **Outputs** : Combined weights cataloged in {subregion}_gene.profile.err.
```
sbatch compute_weights_Amygdala_gene.sh
# and/or
sbatch compute_weights_sACC_gene.sh
```


## 5) Process hg38 GWAS Summary Statistics

Standardizes external GWAS files mapped to the hg38 reference genome. This calculates uniform $Z$-scores across all represented variants, filters out incompatible or missing columns, and generates diagnostic metrics (e.g., $Z$-score and log odds ratio frequency distributions).
```
sbatch run_process-gwas.sh
```


## 6) Apply Weights (Imputation & Association)
Performs expression imputation and calculates trait associations using FUSION.assoc_test.R. Downstream dependencies trigger FUSION.post_process.R to run joint/conditional tests across regions with dense signals, helping distinguish independent eQTL effects from passenger correlations.
```
sh apply_weights.sh
```


## 7) Aggregate Regional TWAS Output
Harvests all calculated association outputs across both the ```sACC``` and ```amygdala``` analytical arms, formatting and consolidating them into an evaluation archive.

- **Outputs** : Comprehensive summary database in ```.Rdata``` format.
```
sbatch run_read_twas_amygdala.sh
# and/or
sbatch run_read_twas_sacc.sh
```


## 8) TWAS plot generation

This script generates a number of different plots and outputs their corresponding tables into `analysis/plots/` and `analysis/tables/`, respectively.

```
sbatch run_generate_twas_plots.sh
```


## References
1. Gusev, A., Ko, A., Shi, H., Bhatia, G., Chung, W., Penninx, B.W., Jansen, R., De Geus, E.J., Boomsma, D.I., Wright, F.A. and Sullivan, P.F., 2016. Integrative approaches for large-scale transcriptome-wide association studies. Nature genetics, 48(3), pp.245-252. https://www.nature.com/articles/ng.3506
2. O’Connell, K.S., Koromina, M., van der Veen, T., Boltz, T., David, F.S., Yang, J.M.K., Lin, K.H., Wang, X., Coleman, J.R., Mitchell, B.L. and McGrouther, C.C., 2025. Genomics yields biological and phenotypic insights into bipolar disorder. Nature, 639(8056), pp.968-975. https://www.nature.com/articles/s41586-024-08468-9  
3. Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118.
