# 🧬 Polygenic Risk Score (PRS) Pipeline

This directory contains the computational workflow for calculating Polygenic Risk Scores (PRS) using **PRSice-2**. This analysis evaluates the shared and distinct genetic liabilities within a cohort of Major Depressive Disorder (MDD) patients, Bipolar Disorder (BP) patients, and healthy controls.

---

## 📊 What is a Polygenic Risk Score?

Unlike single-gene (Mendelian) disorders, psychiatric conditions like MDD and BP are highly polygenic. They are influenced by thousands of single nucleotide polymorphisms (SNPs) across the entire genome, each contributing a small fraction to overall disease risk.

A **Polygenic Risk Score (PRS)** is a single number estimating an individual's genetic liability for a disease. It sums all risk-increasing variants carried by an individual, weighted by their respective effect sizes discovered in massive external genome-wide association studies (GWAS).

Mathematically, the score for individual $i$ is expressed as:

$$PRS_i = \sum_{j=1}^{m} \beta_j \times G_{ij}$$

Where:
* $PRS_i$ is the polygenic score for individual $i$.
* $m$ is the total number of evaluated SNPs.
* $\beta_j$ is the log odds ratio (effect size) of SNP $j$ from the background discovery GWAS summary statistics.
* $G_{ij}$ is the risk allele dosage (0, 1, or 2) carried by individual $i$ at SNP $j$.

---

## 🛠️ Inputs and Outputs

### The Inputs
To execute the scoring pipeline, two distinct genomic structures are required:
1. **Base Data (GWAS Summary Statistics):** External discovery statistics supplying the risk-variant weights (e.g., PGC or UK Biobank datasets).
2. **Target Data (Genotypes):** Internal sample cohort containing the physical DNA calls for the analytical dataset (~750 individuals representing MDD, BP, and healthy control groups), typically formatted as PLINK binary sets (`.bed`, `.bim`, `.fam`).

### The Outputs
Upon completion, the pipeline yields three core components:
1. **Individual Scores:** A sample-by-score matrix detailing the exact polygenic index assigned to every individual in the cohort.
2. **Predictive Fit ($R^2$):** A statistical evaluation detailing how much observed phenotypic variance is accounted for by the genetic score.
3. **P-Value Threshold Matrix:** A resolution table testing multiple statistical stringencies (e.g., $P < 0.05, P < 0.01$) to determine the optimal SNP selection threshold for prediction accuracy.

---

## 🚀 Execution & Command Parameters

The primary execution wrapper is automated through SLURM. It manages core resources and feeds the target data matrices into the scoring engine:

```bash
sbatch 3_run_prsice.sh

```


## Parameter Breakdown:
- --base: Specifies the background public dataset file containing discovery variant weights.

- --target: Points to the local genotype prefix (mdd_bpd_maf01) containing the study cohort.

- --pheno: Path to the phenotype mapping file linking Sample IDs to diagnostic status (e.g., Control, MDD, BP).

- --snp, --chr, --bp, --A1, --A2, --pvalue: Column header mappings matching the internal syntax of the Base GWAS file.

- --stat OR: Directs the software to interpret the input effect size as an Odds Ratio (standard for binary disease models).

- --binary-target T: Explicitly flags the sample outcome matrix as a binary/case-control format rather than a continuous trait.

- --thread 4: Allocates multithreading processing loops, explicitly locked to match the allocated cluster core headers (#SBATCH -c 4).