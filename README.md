# High-Resolution Transcriptomic Analysis of Bipolar Disorder vs. Major Depressive Disorder

A computational pipeline investigating the distinct molecular signatures of **Bipolar Disorder (BP)** and **Major Depressive Disorder (MDD)** within the human **Amygdala** and **subgenual Anterior Cingulate Cortex (sgACC)** across multiple resolution levels: Gene, Exon, Junction, and Transcript (Isoform).

---

## 📌 Project Overview & Rationale

### Why BP vs. MDD?
Distinguishing between MDD and BP remains a significant clinical challenge due to overlapping depressive symptoms. Identifying distinct genomic boundaries is critically important to:
* **Improve Early Differential Diagnosis:** Reducing misdiagnosis rates and lag time to correct treatment [[1, 2](#references)].  
* **Guide Personalized Treatment:** Avoiding the risk of antidepressant-induced mania in BP patients [[3, 4](#references)].
* **Map Disease Biology:** Uncovering the unique biological mechanisms underlying both conditions [[5](#references)].

## 🧠 Brain Region Breakdown

### Amygdala: The Emotional Sentinel
The amygdala acts as a critical hub for determining the emotional significance of stimuli (salience) and forming emotional memories.
* **Role in MDD:** Demonstrates **hyperactivity** and heightened reactivity to negative emotional stimuli, coupled with a failure of top-down regulation from the prefrontal cortex (PFC).
* **Role in BP:** Functions in a **state-dependent** manner. It shows profound hyperactivity/increased volume during mania, but variable or reduced activity during depressive phases.

### sgACC: The Mood Regulator
The sgACC (Brodmann Area 25) modulates emotional, visceral, and neuroendocrine responses (such as cortisol regulation).
* **Shared Structural Abnormality:** Both MDD and BP patients show a consistent reduction in gray matter volume in the sgACC, typically driven by a loss of glial cells rather than neurons.
* **Metabolic Hyperactivity:** During active depressive episodes in both disorders, the sgACC displays elevated metabolic activity ("hot spot") that acts as a therapeutic target for Deep Brain Stimulation (DBS).

### Functional Profiles at a Glance

| Feature | Major Depressive Disorder (MDD) | Bipolar Disorder (BP) |
| :--- | :--- | :--- |
| **Amygdala Activity** | Consistently **Hyperactive** (especially to negative stimuli). | **State-Dependent** (Hyperactive in mania, variable/hypoactive in depression). |
| **sgACC Structure** | **Reduced Volume** (Trait-like abnormality). | **Reduced Volume** (Trait-like abnormality, similar to MDD). |
| **sgACC Activity** | **Hyperactive** during the depressive state. | **Hyperactive** during the depressive state. |
| **Connectivity** | Disrupted connectivity between the amygdala and frontal regulatory areas; increased adolescent sgACC-amygdala connectivity. | Distinct patterns of **prefrontal-amygdala** and **sgACC-amygdala** connectivity, particularly during pronounced regulatory failure in mania. |

---

## 🔍 Why Resolution Matters: Transcriptomics Beyond the Gene

Standard gene-level analysis often masks crucial distinctions because MDD and BP share risk genes. Investigating finer-grained features—**Exons, Junctions, and Transcripts (Isoforms)**—allows us to bypass these limitations through the lens of **Alternative Splicing**.

* **Differential Transcript (Isoform) Usage:** A single gene can produce multiple protein versions with opposing functions. MDD and BP may express the same gene, but produce entirely different isoform profiles.
* **Exon & Junction Level Analysis:** RNA processing and splicing alterations are fundamental molecular risks in neuropsychiatric conditions. Mapping splice junctions catches specific regulatory defects that total gene expression metrics completely miss.
* **Distinct Biological Pathways:** High-resolution sequencing reveals that while BP pathways are dominated by **immune system activation and immune-regulatory factors**, MDD pathways lean heavily toward **stress response and metabolic dysregulation**.

---

## 🛠️ Project Pipeline & Methods

The structured analytical workflow implemented in this repository is outlined below. 

Click the links below to access the code, documentation, and specific results for each step of the pipeline:

* [**Differential Expression Analysis:**](/Diffierentially_expression_feature) Identification of Differentially Expressed Features (**DEFs**) across genes, transcripts, exons, and junctions.
* [**WGCNA Network Analysis:**](/WGCNA) Construction of weighted gene co-expression networks to find highly correlated feature modules.
* [**Gene-Set Enrichment:**](/Gene_set_enrichmnet) Functional enrichment and pathway mapping of identified DEFs.
* [**Fusion-TWAS:**](TWAS) Transcriptome-Wide Association Studies integrating GWAS summary statistics with reference expression weights.
* [**IsoTWAS:**](isoTWAS) Isoform-specific TWAS evaluations focusing on isoform-specific expression (ISE) and tissue-specific regulation.
* [**Leafcutter Splicing Analysis:**](Leafcutter) Annotation-free quantification of RNA splicing variation focusing explicitly on intron usage across samples.

---

## 📚 References

1. Panagiotaropoulou, G. et al. (2025) ‘Identifying genetic differences between bipolar disorder and major depression through multiple genome-wide association analyses’, *The British Journal of Psychiatry*, 226(2), pp. 79–90. doi:10.1192/bjp.2024.125.
2. Yang, R. et al. (2023) ‘Differentiation between bipolar disorder and major depressive disorder in adolescents: from clinical to biological biomarkers’, *Frontiers in Human Neuroscience*, 17, p.1192544.
3. Tonozzi, T. R. et al. (2018) ‘Pharmacogenetic Profile and Major Depressive and/or Bipolar Disorder Treatment: a Retrospective, Cross-Sectional Study’, *Pharmacogenomics*, 19(15), pp. 1169–1179. doi:10.2217/pgs-2018-0088.
4. Wilcox, C. (2020) ‘Can Genetics Help Us Distinguish Bipolar Disorder from Major Depression?’, *NEJM Journal Watch*, p.NA51036.
5. Panagiotaropoulou, G. et al. (2025) ‘Identifying genetic differences between bipolar disorder and major depression through multiple genome-wide association analyses’, *The British Journal of Psychiatry*, 226(2), pp.79-90.