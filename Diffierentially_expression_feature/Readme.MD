# Differential Expression Feature (DEF) Pipeline via qSVs

This directory contains the computational workflow for evaluating Differential Expression Features (DEFs) between Major Depressive Disorder (MDD) and Bipolar Disorder (BP) cohorts. To account for biological heterogeneity and technical artifacts inherent to postmortem human brain tissue, this framework utilizes **Quality Surrogate Variables (qSVs)** alongside conventional covariates.

The pipeline spans four structural genomic resolution layers: **Genes, Exons, Splice Junctions, and Transcripts (Isoforms)**, processing the **amygdala** and **subgenual anterior cingulate cortex (sgACC)** architectures independently.

---

## 🛠️ Pipeline Execution & Script Breakdown

### 1) Orchestrate DE Pipeline via SLURM (`qSV_DE_MDD_BP.sh`)

Serves as the high-throughput cluster wrapper responsible for allocating computing resources ($50\text{ GB}$ RAM) and initiating the full differential expression engine. It standardizes environmental configurations via Conda modules to ensure exact statistical reproducibility across executions.

* **Primary Function:** Invokes `qSV_DE_MDD_BP.R` within a dedicated node allocation.



```bash
# Submit the complete multi-omic DE execution array:
sbatch qSV_DE_MDD_BP.sh

```

---

### 2) Core Multi-Resolution Differential Modeling (`qSV_DE_MDD_BP.R`)

Acts as the central engine of the differential expression framework. The script normalizes transcriptomic profiles, handles covariate matrices, and executes parallelized regression modules across tissues and feature resolutions. It sequentially completes the following tasks:

* **Feature Standardization:** Loads localized RangedSummarizedExperiment (`rse`) structures and standardizes metadata definitions across layers by enforcing absolute naming conventions (`common_feature_id`, `common_gene_symbol`, and `common_gene_id`).


* **Covariate Layering:** Pairs anatomical subsets with explicit mathematical designs. Models isolate primary diagnoses (`PrimaryDx`) while dynamically regressing out demographic features (`Age`, `Sex`), technical matrix statistics (`mitoRate`, `rRNA_rate`, `RIN`), population ancestry dimensions (`snpPCs`), and tissue degradation layers derived from quality surrogate variable matrices (`qSV_mat`).


* **High-Throughput Testing Matrix:** Leverages custom wrapper protocols to compute linear fits across all four distinct levels of biological resolution:


* **Gene-Level Resolution:** Evaluates aggregate count variances mapped across full gene structures.


* **Exon-Level Resolution:** Screens for regional signal imbalances to signal potential exon skipping or alternative inclusion.


* **Junction-Level Resolution:** Quantifies absolute count boundaries at splicing coordinates to record active splice variations.


* **Transcript-Level Resolution (Isoform):** Investigates precise full-length transcript usage matrices evaluated via TPM log transformations ( $log_2(\text{TPM} + 1)$ ) rather than raw counts.




* **Outputs:** Consolidates multi-omic statistical tables into compressed summary archives (`BP_MDD_All_results.RDS`) containing $F$-statistics, localized $P$-values, and False Discovery Rate (FDR) $q$-values for both target subregions.



---

### 3) Core Splicing & Linear Regression Backends (`run_DE.R`)

Supplies the underlying mathematical frameworks used by the master script to calculate feature variances. It leverages the robust `limma` and `edgeR` ecosystems to process variance trends:

* `get_eBayes()`: Dynamically branches according to the input configuration. For counts-based resolutions (Gene, Exon, Junction), it calculates scaling elements via `edgeR::calcNormFactors` and structures empirical weights via `limma::voom` to control for the inherent mean-variance relationship of RNA-Seq data. For transcript configurations, it bypasses voom modeling, directly implementing linear fits on log-transformed expression targets. It then evaluates final models via empirical Bayes moderation (`limma::eBayes`).


* `run_DE()`: Extracts individual design coefficient comparisons across diagnostic boundaries. Because the model accommodates complex, multi-component diagnostic profiles, it evaluates robust $F$-statistics alongside standard $t$-statistics, applies Benjamini-Hochberg FDR adjustments to control multiple testing error rates, and compiles scannable significance tables.



```
# Internal helper package called automatically inside qSV_DE_MDD_BP.R
source("codes/run_DE.R")

```