# Gene-Set Enrichment Analysis (GSEA) & Pathway Mapping Pipeline

This directory contains the computational workflow for conducting functional gene-set enrichment testing on the Differential Expression Features (DEFs) identified between Major Depressive Disorder (MDD) and Bipolar Disorder (BP) cohorts. By mapping statistically significant features across multiple structural resolutions (Genes, Exons, Junctions, and Transcripts) and directional splits (Up-regulated, Down-regulated, and Whole list), this pipeline uncovers the underlying biological pathways and molecular functions operating within the **amygdala** and **subgenual anterior cingulate cortex (sgACC)**.

---

## 🛠️ Pipeline Execution & Script Breakdown

### 1) Orchestrate Enrichment Analysis via SLURM (`Enrichment.sh`)

Acts as the high-performance cluster wrapper script designed to allocate $20\text{ GB}$ of RAM and launch the core enrichment routines. It configures the runtime environment using standard Conda modules to ensure exact statistical reproducibility across execution iterations.

* **Primary Function:** Spawns `Enrichment.R` inside a single-core cluster allocation.

```bash
# Submit the complete functional pathway profiling pipeline:
sbatch Enrichment.sh

```

---

### 2) Core GO Ontology Profiling Engine (`Enrichment.R`)

Serves as the main analytical pipeline for translating multi-resolution differential expression outputs into localized biological profiles. The script utilizes a standardized background universe and handles data translation, partitioning, and automated plotting. It executes the following steps:

* **Background Universe Construction:** Loads global differential results (`BP_MDD_All_results.RDS`), compiles all tested feature IDs, maps them across ENSEMBL categories, and utilizes the `org.Hs.eg.db` database to resolve them into a standardized background of Entrez gene identifiers.
* **Significance Sorting & Stratification:** Employs custom extraction modules (`get_signif`) to subset target features by subregion based on an FDR adjusted threshold ($q < 0.05$). It isolates directional signals into three independent profiles:
* **Whole List Matrix:** Maps all statistically significant features regardless of fold-change direction.
* **Up-Regulated Cohort:** Captures genes and isoforms demonstrating enhanced expression in Bipolar Disorder relative to MDD.
* **Down-Regulated Cohort:** Identifies features demonstrating suppressed expression in Bipolar Disorder relative to MDD.


* **Hypergeometric Enrichment Testing:** Flattens tissue and feature intersections to avoid redundancy and runs over-representation analysis using `clusterProfiler::enrichGO`. It simultaneously models across all primary Gene Ontology (GO) frameworks—**Biological Process (BP)**, **Cellular Component (CC)**, and **Molecular Function (MF)**—using Benjamini-Hochberg multi-test adjustments.
* **Automated Visualization:** Automatically filters out null gene sets and generates detailed pathway dotplots (`clusterProfiler::dotplot`) displaying the top 25 enriched terms.
* **Outputs:** Saves separate raw results data files (`enrichGO_whole_list.rda`, `enrichGO_up.rda`, `enrichGO_down.rda`) alongside high-resolution diagnostic dotplot images (`.jpg`) and combined multi-page reference documents (`.pdf`) under the project graphs tree.