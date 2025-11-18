# Comparative Transcriptomics of Major Depressive Disorder (MDD) and Bipolar Disorder (BP) in the Amygdala and sACC

## The Challenge of Differential Diagnosis in Affective Disorders

Major Depressive Disorder (MDD) and Bipolar Disorder (BP), particularly during episodes of Bipolar Depression (BPD), present a significant challenge in psychiatric differential diagnosis due to extensive symptomatic overlap. This ambiguity often leads to high rates of misdiagnosis, with estimates suggesting that up to 40% of patients with BPD are initially diagnosed with MDD  [[1](#references)]. Such diagnostic error is clinically catastrophic; conventional antidepressant monotherapy, often prescribed for MDD, can exacerbate the course of BP, potentially triggering rapid cycling or manic episodes. Consequently, the lack of objective, validated molecular biomarkers hinders the ability to stratify patients accurately, delaying the implementation of appropriate mood stabilizers or antipsychotic regimens required for BP management [[2](#references)].   

The reliance solely on behavioral and self-reported symptoms underscores the pressing need for objective molecular stratification. Large-scale comparative transcriptomic analyses serve to move beyond descriptive symptoms by systematically comparing the gene expression profiles of MDD and BP patients. The overarching goal is to identify reproducible biological distinctions that can inform prognosis, predict therapeutic response, and ultimately validate novel targets for mechanism-based pharmacological intervention [[1](#references)] . Establishing these molecular disjunctions is foundational to mitigating the high rates of morbidity, disability, and suicide risk associated with these chronic affective disorders[[2](#references)].   

## The Strategic Significance of Amygdala and sACC

The molecular investigation of psychiatric disorders demands a focus on brain regions central to the disease pathophysiology. The Amygdala and the subgenual Anterior Cingulate Cortex (sACC) are strategically chosen for comparative transcriptomic analyses because they constitute the core nodes of the limbic-cortical regulatory circuit that governs mood, emotion processing, and valuation [[3](#references)] . Functional and structural neuroimaging studies have consistently implicated both regions in the pathophysiology of MDD and BP, confirming their status as primary pathological hubs. The Amygdala, for instance, plays a crucial role in threat detection and emotional salience, with its enhanced activity and altered connectivity widely associated with stress-induced anxiety and depressive behaviors [[4](#references)].   

Conversely, the sACC is critical for the executive regulation of mood, and transcriptional changes within this region are considered highly relevant to MDD and BP pathology [[2](#references)]. Transcriptomic analysis in these specific post-mortem brain regions—as opposed to heterogeneous cortical samples or peripheral blood—provides the highest biological fidelity by capturing the molecular consequences directly linked to the observed functional and structural dysregulation of mood circuitry. Studies focusing on these regions, such as a large-scale analysis involving 846 samples across 458 individuals in the sACC and Amygdala, demonstrate that these specific anatomical locations are indispensable for providing the initial mechanistic understanding of these disorders and highlighting potential drug targets[[5](#references)].   

## Methods used in this project

![Flowchart](https://github.com/karbalaei/MDD_BP/blob/main/graphs/Flowchart.jpg) 

In this proect , these analysis has been done yet :

- A- Find differentially expressed features( gene, transcripts, exon, junction) or **DEFs**.
- B- Co-expressiopn networks analyzed using **WGCNA** package
- C- **Gene-set enrichmnt** of DEFs.
- D- Functional Summary-data ImputatiOn Transcriptome-Wide Association Studies (**Fusion-TWAS**) which is primarily designed to use GWAS summary statistics combined with expression weights calculated from a separate reference panel of individuals with both genotype and expression data. 
- E- **IsoTWAS** which goes beyond gene-level expression and focuses on isoform-specific expression (ISE). This is crucial because different isoforms of the same gene can have distinct functions and tissue-specific regulation.
- F- **Leafcutter** which is a computational tool designed to analyze splicing variation across samples and tissues, specifically focusing on intron usage.

## References
1. Malik, S., Singh, R., Arora, G., Dangol, A. and Goyal, S., 2021. Biomarkers of major depressive disorder: knowing is half the battle. Clinical Psychopharmacology and Neuroscience, 19(1), p.12.
2. Goes, F.S., Collado-Torres, L., Zandi, P.P., Huuki-Myers, L., Tao, R., Jaffe, A.E., Pertea, G., Shin, J.H., Weinberger, D.R., Kleinman, J.E. and Hyde, T.M., 2025. Large-scale transcriptomic analyses of major depressive disorder reveal convergent dysregulation of synaptic pathways in excitatory neurons. Nature communications, 16(1), p.3981.
3. Hu, P., Lu, Y., Pan, B.X. and Zhang, W.H., 2022. New insights into the pivotal role of the amygdala in inflammation-related depression and anxiety disorder. International journal of molecular sciences, 23(19), p.11076.
4. Ramasubbu, R., Konduru, N., Cortese, F., Bray, S., Gaxiola-Valdez, I. and Goodyear, B., 2014. Reduced intrinsic connectivity of amygdala in adults with major depressive disorder. Frontiers in psychiatry, 5, p.17.
5. Jaffe, A.E., Tao, R., Page, S.C., Maynard, K.R., Pattie, E.A., Nguyen, C.V., Deep-Soboslay, A., Bharadwaj, R., Young, K.A., Friedman, M.J. and Williamson, D.E., 2022. Decoding shared versus divergent transcriptomic signatures across cortico-amygdala circuitry in PTSD and depressive disorders. American Journal of Psychiatry, 179(9), pp.673-686.