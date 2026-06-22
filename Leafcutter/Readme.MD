# Leafcutter Alternative Splicing Analysis Pipeline 

This directory contains the computational pipeline for conducting annotation-free **Differential Splicing (DS)** analysis between Major Depressive Disorder (MDD) and Bipolar Disorder (BP) cohorts using Leafcutter. Instead of relying on full-length transcript reconstructions, Leafcutter estimates alternative splicing variations by quantifying **intron excision clusters** directly from RNA-Seq split-read mappings across the **amygdala** and **subgenual anterior cingulate cortex (sgACC)**.

---

## 🛠️ Pipeline Execution & Steps

### Step 1: Model Construction & Sample Alignment

Prepares the workspace and parses the experimental design metadata for Leafcutter.

* **R Component (`leafcutter_step_1_copy_files_create_models.R`):** Loads the junction metadata file (`rse_jxn.Rdata`), filters it to isolate MDD and Bipolar samples, and formats the diagnostic comparisons. It integrates sample degradation covariates using Quality Surrogate Variables (`qSV_mat.Rdata`) alongside demographic factors (Age, Sex), technical metrics (RIN, alignment rates), and population structure vectors (`snpPCs`). It then exports plain, headerless covariate maps for both the Amygdala and sACC tissues as required by Leafcutter's command-line interface. Finally, it handles file system operations to stage raw junction `.bed` fragments.
* **SLURM Automation Wrapper (`leafcutter_step_1_copy_files_create_models.sh`):** Executes the baseline R script and builds a centralized text index mapping the absolute file system address of every valid junction target.
* **Outputs:** Regional design models (`modSep_Amygdala_for_leafcutter.txt`, `modSep_sACC_for_leafcutter.txt`) and a primary registry index file (`leafcutter_DS_juncfiles_covariates.txt`).

```bash
# Execute dataset collection and structural metadata modeling:
sbatch leafcutter_step_1_copy_files_create_models.sh

```

---

### Step 2: High-Performance Intron Clustering

Executes the intron clustering stage to map overlapping split-reads across the unified sample registry.

* **SLURM Wrapper (`leafcutter_step_2_Differential_Splicing_BP_vs_MDD_clustering.sh`):** Requests a high-memory cluster allocation ($100\text{ GB}$ RAM) to pool reads and process intermediate clustering operations across all samples simultaneously.
* **Clustering Backend Core (`leafcutter_cluster_regtools.py`):** An optimized execution architecture designed to parse millions of split-read records. It scans the raw junction files compiled in Step 1, isolates high-quality biallelic mappings, and drops ambiguous read trends or extreme structural outliers via strict thresholds ($30$ minimum reads per cluster, maximum intron bridge span of $500,000\text{ bp}$). Overlapping intervals are grouped using graph-based link-refining routines to construct distinct, sample-by-junction count grids.
* **Outputs:** Compressed multi-column matrices tracking read mappings across junctions and individuals (`leafcutter_perind_numers.counts.gz`).

```bash
# Generate unified splice junction clusters:
sbatch leafcutter_step_2_Differential_Splicing_BP_vs_MDD_clustering.sh

```

---

### Step 3: Differential Intron Excision Testing

Tests for differential splicing between MDD and BP phenotypes within the localized tissue targets.

* **Subregion Execution Workers (`leafcutter_step_3_DS_Amygdala.sh` and `leafcutter_step_3_DS_sACC.sh`):** Load the Leafcutter application layer and run the core differential analysis code (`leafcutter_ds.R`) using $4$ threads. The scripts combine the global sample count matrix generated in Step 2 with the specific covariate text configurations produced in Step 1. This uses a Dirichlet-multinomial generalized linear model to test for changes in intron excision ratios across diagnostic conditions.
* **Master Controller Wrapper (`leafcutter_step_3.sh`):** A bash script used to concurrently submit both tissue evaluation scripts to the cluster queue.
* **Outputs:** Tissue-specific statistical results archives containing effect sizes ($\Delta\text{PSI}$ values) and multi-test adjusted $P$-values saved to `Amygdala_covariates` and `sACC_covariates` prefixes.

```bash
# Launch concurrent differential splicing tests across both subregions:
sh leafcutter_step_3.sh

```

---

### Step 4: Download & Build Transcriptome Annotation Maps

Prepares reference layers to map and translate raw splicing coordinates into biological features.

* **SLURM Script (`leafcutter_step_4_download_hg38_prepare_annotation.sh`):** Downloads the human genome reference database (GENCODE v25 annotation file mapped to hg38 coordinates). It modifies user access permissions for downstream scripts and launches a Perl extraction routine (`gtf2leafcutter.pl`). This utility translates the global structural database into localized intron, exon, and gene coordinate reference tracks.
* **Outputs:** Standardized exon-to-gene mapping profiles and transcript boundary references (`gencode_hg38_all_introns.bed.gz`, `gencode_hg38_threeprime.bed.gz`, etc.) under the project results workspace.

```bash
# Fetch hg38 annotations and build reference maps for results conversion:
sbatch leafcutter_step_4_download_hg38_prepare_annotation.sh

```

### Step 5: Post-Testing Annotation & Data Packaging

Converts raw intron genomic coordinates into biologically interpretable datasets by cross-referencing external databases.

* **SLURM Execution Matrix (`leafcutter_step_5_Amygdala.sh` & `leafcutter_step_5_sACC.sh`):** Allocates cluster nodes ($20\text{ GB}$ RAM) to handle data compilation. The scripts execute Leafcutter's post-processing utility (`prepare_results.R`) independently for each subregion.


* **Annotation Engine Backend (`prepare_results.R`):** Merges raw cluster significance datasets and junction effect-size counts (`leafcutter_perind_numers.counts.gz`). Using an adjusted significance threshold ($\text{FDR} < 0.05$), it queries the `gencode_hg38` reference tracks constructed in Step 4 to determine splicing classifications. It catalogs junctions as **annotated** or maps novel variations into structural categories: *cryptic five-prime, cryptic three-prime, cryptic unanchored*, or *novel annotated pairs*. It also performs Principal Component Analysis (PCA) on intron excision ratios.


* **Outputs:** Consolidated data structures (`step_5_Amygdala.Rdata` and `step_5_sACC.Rdata`) storing clean sample metadata, junction ratio mappings, and PCA summaries.



```bash
# Submit annotation formatting configurations for each tissue arm:
sbatch leafcutter_step_5_Amygdala.sh
sbatch leafcutter_step_5_sACC.sh

```

---

### Step 6: Multi-Omic Filtering & Structural Plotting

Integrates differential expression results with alternative splicing signals to generate localized visual transcripts and tables.

* **SLURM Automators (`leafcutter_step_6_Amygdala.sh` & `leafcutter_step_6_sACC.sh`):** Handle execution environment settings via Conda dependencies to trigger regional downstream scripting modules.


* **Visual Synthesis Engines (`leafcutter_step_6_Amygdala.R` & `leafcutter_step_6_sACC.R`):** Load the RData packages compiled in Step 5 alongside global differential expression tables (`BP_MDD_All_results.RDS`). The scripts dynamically isolate significant Differentially Expressed Genes (DEGs; $\text{FDR} < 0.05$) and tag overlapping intron clusters. The workflow iterates through all significant features using nested `purrr::walk` mapping loops to generate outputs across three graphical domains:


1. **Cluster Count Summaries:** Generates detailed junction statistics, tracking percent spliced-in values ($\Delta\text{PSI}$).


2. **Intron Excision Cluster Plots:** Renders smooth splice-graphs where curve paths match the mean mapping percentages for each clinical group.


3. **Gene-Model Transcripts:** Renders mapped splicing clusters against standard black box-and-whisker exon boundaries.




* **Outputs:** Exports clean final spreadsheets (`Amygdala_sig_clusters.csv`, `Amygdala_sig_introns.csv`) and structured vector graphic transcripts (`_graph.pdf` and `_gene.pdf`) into the project graphs tree.



```bash
# Execute visual rendering workflows for calculated features:
sbatch leafcutter_step_6_Amygdala.sh
sbatch leafcutter_step_6_sACC.sh

```

---

### Shared Plotting Backends (`make_*.R`)

These background scripts are sourced directly by Step 6 to construct tables and graphical layouts:

* `make_cluster_plot_Final.R`: Uses `ggplot2::geom_curve` and custom scale palettes to overlay split-read arc trajectories. It separates annotated reference tracks from cryptic structural boundaries using variable line thicknesses.


* `make_gene_plot_Final2.R` / `make_gene_plot_Final.R`: Matches variable intron positions against absolute chromosomal windows, handling coordinate resizing to ensure smooth visualization of compact elements alongside massive non-coding fields.


* `make_cluster_table_Final.R`: Constructs minimal matrix frames via `gridExtra::tableGrob` and handles layout alignment adjustments to generate clean tabular keys.