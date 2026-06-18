# 🧬 What is isoTWAS?

- **The Concept**: Traditional TWAS evaluates associations at the aggregate gene level, which can obscure critical disease mechanisms. **isoTWAS** extends this framework by shifting the resolution to **isoform-specific expression (ISE)** and transcript-level features[[1](#references)]. This approach is uniquely powerful for distinguishing complex psychiatric conditions like MDD and BP, where alternative splicing and variable transcript usage often drive pathobiology even when total gene-level expression remains unchanged[[2](#references)].

- **Prediction Models** : Similar to gene-level frameworks, isoTWAS trains predictive models (such as Elastic Net, LASSO, and Blended architectures) to estimate the genetically regulated component of expression—but does so independently for individual transcripts or isoforms[[1](#references)].

- **Association Step**: By integrating transcript-specific predictive weights ($\mathbf{w}$) with major GWAS summary statistics (e.g., PGC or 23andMe) and an LD reference panel, isoTWAS computes isoform-level association $Z$-scores[[1](#references)]. This identifies specific transcript variations that are significantly correlated with phenotypic risk, allowing for high-resolution molecular mapping[[1 , 2 ](#references)].

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

## 2) Build Isoform Expression Matrices & Target Gene Catalogs (`2_build_expression_matrix.R`)

Constructs the core transcript expression architectures required for downstream modeling. The script dynamically imports the localized sample expression data (`rse_tx`) alongside pre-calculated Transcript PCs. Confounding variations—including diagnostic status (`PrimaryDx`), genetic ancestry (`snpPCs`), sex, and transcript-level batch effects—are rigorously controlled using a linear framework via `jaffelab::cleaningY` to isolate **"clean expression" residuals**.

To maximize analytical power, the workflow loads prior gene-level heritability files from the companion TWAS analysis, filtering out feature windows that lack statistical significance ($P > 0.05$ or $h^2 \le 0$). It isolates genes possessing **multiple distinct transcript isoforms** ($n > 1$) and maps them against the significance index.

- **Outputs:** A clean, standardized transcript-level expression matrix saved to `expression_matrix.rda` alongside subregion-specific batch maps (`gene_list_isotwas_{region}.txt`) containing chromosome coordinates and gene indexes for distributed cluster scheduling.

```bash
# Run to generate matrices:
Rscript 2_build_expression_matrix.R --region Amygdala
# and/or
Rscript 2_build_expression_matrix.R --region sACC

```

---

### 3) Compute Distributed Isoform Weights via Distributed Array Jobs (`3_compute_isotwas_*.sh`)

Automates high-throughput parallel computation of isoform-specific genetic prediction models across thousands of features using the SLURM workload manager. The shell wrappers read the target gene catalogs compiled in Step 2, parse the row indexes corresponding to the active `SLURM_ARRAY_TASK_ID`, and spawn dedicated `Rscript Compute_isotwas.R` threads. Each task dynamically invokes predictive regularized modeling algorithms to build custom genetic weights for every transcript variant within the locus.

* **Outputs:** Independent, localized transcript-level predictive weight models mapped across $6,499$ arrays for the Amygdala and $7,100$ arrays for the sACC.

```bash
# Submit highly parallelized job arrays to the cluster:
sbatch 3_compute_isotwas_Amygdala.sh
# and/or
sbatch 3_compute_isotwas_sACC.sh

```

---

### 4) Construct Reference Linkage Disequilibrium Matrices (`4_create_LD_matrix.R`)

Generates structural linkage disequilibrium (LD) covariance references across the genome using raw European reference panels from the 1000 Genomes Project (`1000G.EUR.*.bed`). The pipeline systematically iterates through all autosomes (Chromosomes 1–22), mapping the genotype markers using `bigsnpr::snp_readBed2` and evaluating absolute correlation matrices via `bigsnpr::snp_cor`. Variant identifiers are synchronized across rows and columns using unique marker IDs to allow seamless lookup during the final association steps.

* **Outputs:** 22 autosome-specific covariance matrices saved to `isotwas/LD_matrix/chromosome_{1:22}.rds`, supplying the necessary linkage framework to compute isoform-level association $Z$-scores with external GWAS summary statistics.

```bash
# Build background LD covariance layers:
Rscript 4_create_LD_matrix.R

```

---

### 5) Format & Process hg38 GWAS Summary Statistics (`5_process_gwas_isotwas.R`)

Standardizes and cleans external Psychiatric Genomics Consortium (PGC) and 23andMe major depression summary statistics mapped to the hg38 genome reference block. The script enforces a rigorous column layout, re-indexes variant identifiers against the unique project `.bim` coordinates generated in Step 1, and back-calculates missing statistical parameters. Specifically, it derives uniform alignment $Z$-scores from log odds ratios and standard errors ($Z = \frac{\text{LogOR}}{\text{StdErrLogOR}}$) and estimates effective sample sizes ($N$) dynamically per variant using minor allele frequencies.

* **Outputs:** A highly optimized, chromosome-annotated genomic file saved as `PGC_UKB_23andMe_depression_genome-wide_isotwas_hg38_ver2.txt` ready for cross-referencing.

```bash
# Process and standardize the trait summary statistics:
Rscript 5_process_gwas_isotwas.R

```

---

### 6) Compute Isoform Association Tests via Burden Modeling (`6_running_burden_test.R`)

Carries out the core summary-based association testing across all 22 autosomes using the `isotwas::burdenTest` framework. The script matches the pre-computed, transcript-specific prediction weights generated in Step 3 against the processed GWAS dataset (Step 5) and the chromosome-specific reference linkage disequilibrium matrices (Step 4). For features meeting the minimum predictive threshold ($R^2 > 0.01$), it calculates a weighted burden $Z$-score and an association $P$-value, capturing whether genetically predicted alternative transcript variations drive disease liability. For highly significant loci, it triggers localized permutation testing ($1,000$ iterations) to provide robust empirical $P$-values.

* **Outputs:** Comprehensive subregion burden tables cataloging transcript-level association profiles across the genome, saved as `_BurdenTest_results_MDD.RDS` and `.txt` files.

```bash
# Launch the burden testing arrays across the cluster:
sbatch 6_running_burden_test_Amygdala.sh
# and/or
sbatch 6_running_burden_test_sACC.sh

```

---

### 7) Multi-Stage Multiple Testing Correction & FDR Filtering (`7_FDR_calculaton_isotwas.R`)

Applies a specialized hierarchical multiple-testing correction framework designed specifically for multi-transcript isoform settings. To control the False Discovery Rate (FDR) without sacrificing power due to tight transcript co-regulation, the pipeline performs a **two-stage aggregation and confirmation test**:

1. **Screening Phase:** Evaluates aggregate gene-level risk by pooling individual transcript $P$-values via the Cauchy combination test (`isotwas::p_screen`) and applies Benjamini-Hochberg FDR adjustment ($\alpha_1 = 0.05$) to find significant overall loci.
2. **Confirmation Phase:** For genes passing the global screening tier, individual isoforms are evaluated through a localized conditional confirmation threshold (`isotwas::p_confirm`) to pin down the exact causal transcripts.

Finally, the script cross-references genomic bounds to isolate significant elements located within $1\text{ Mb}$ windows that may warrant downstream fine-mapping.

* **Outputs:** Final prioritized target tables saved as `_all_results_MDD.csv` and filtered, high-confidence files containing biologically meaningful transcripts (`_FDR_filtered_results_MDD.csv`).

```bash
# Run multiple-testing corrections across the calculated features:
Rscript 7_FDR_calculaton_isotwas.R --region "Amygdala" --diag "MDD"
# and/or
Rscript 7_FDR_calculaton_isotwas.R --region "sACC" --diag "MDD"

```

## References
1- Bhattacharya, A. et al. (2022) 'isoTWAS: Transcript-level transcriptome-wide association studies via structured lasso and blended learning', Genetics, 222(4), p.iyac145. doi:10.1093/genetics/iyac145.

2- Gandal, M. J. et al. (2018) 'Transcriptome-wide isoforms-level imbalances in the brains of individuals with severe mental illnesses', Science, 362(6420), p.eaat8127. doi:10.1126/science.aat8127.