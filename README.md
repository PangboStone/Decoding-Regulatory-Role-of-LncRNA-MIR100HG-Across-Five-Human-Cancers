# Decoding Regulatory Role of LncRNA MIR100HG Across Five Common Human Cancers
## Abstract
This project investigated MIR100HG-regulated mechanisms across five common cancers by integrating patient stratification, gene expression, methylation, and transcription factor–target relationships.

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

## 1. Executive Summary
This project presents a comprehensive computational biology workflow to investigate the regulatory landscape of the long non-coding RNA (lncRNA) MIR100HG across five distinct cancer types: Lung Adenocarcinoma (LUAD), Pancreatic Adenocarcinoma (PAAD), Prostate Adenocarcinoma (PRAD), Skin Cutaneous Melanoma (SKCM), and Stomach Adenocarcinoma (STAD). By integrating multi-omics data including gene expression (RNA-seq) and DNA methylation, the analysis explores two primary regulatory axes. Firstly, it examines the potential auto-regulation of MIR100HG through the methylation status of its own promoter region. Secondly, it delineates the downstream regulatory network by identifying key target genes and their interactions. To refine this network, a sophisticated ensemble model of Graph Attention Networks (GAT) is developed to predict and validate high-confidence regulatory relationships, moving beyond simple correlation to a more robust, machine learning-driven conclusion.

## 2. Problem Statement
Long non-coding RNAs (lncRNAs) are emerging as critical regulators in various cellular processes, and their dysregulation is frequently implicated in cancer pathogenesis. MIR100HG, the host gene for the miR-100/let-7a-2/miR-125b-1 cluster, has been identified as a significant player in tumorigenesis, yet its precise regulatory mechanisms and functional network remain incompletely understood. This project aims to bridge this gap by systematically dissecting the regulatory functions of MIR100HG. The central objectives are:
- To investigate the relationship between MIR100HG's promoter DNA methylation and its own expression, exploring a potential epigenetic self-regulation loop.
- To identify and characterize the downstream regulatory sub-network controlled by MIR100HG.
- To construct a high-confidence MIR100HG-centric regulatory network by leveraging advanced graph neural network models to reduce noise and infer robust interactions.

## 3. Tech Stack & Environment
- **Language:** Python 3.1
- **Core Libraries:** Pandas, NumPy, SciPy
- **Data Analysis & Statistics:** scikit-learn, statsmodels
- **Network Biology:** NetworkX
- **Visualisation:** Matplotlib, Seaborn
- **Machine Learning:** PyTorch, PyTorch Geometric
- **Graph Embeddings:** node2vec

## 4. Dataset
- All research dataset was sourced from The Cancer Genome Atlas (TCGA).
- Genomic annotations and probe maps were obtained from UCSC Genome Browser resources.
This study utilises publicly available multi-omics data, primarily sourced from The Cancer Genome Atlas (TCGA) repository. The datasets for each of the five cancer types (LUAD, PAAD, PRAD, SKCM, STAD) include:
### Data Types
- **Cancer Gene Expression (TCGA):** Normalised gene expression levels ($log2(TPM+0.001)$) for five cancer types: PAAD, LUAD, SKCM, PRAD, and STAD.
- **Normal Tissue Gene Expression (GTEx):** Normalised gene expression levels for five corresponding normal tissues: pancreas, lung, skin, prostate, and stomach.
- **Transcription Factor-Target Associations (ENCODE):** A curated list of associations between TFs and their target genes.
- **Clinical and Survival Data (TCGA Pan-Cancer):** Patient clinical information, including histological subtype and survival endpoints (OS, DSS, PFI).
- **Gene Methylation Data (TCGA, 450k):** Genome-wide methylation values for the five cancer types. *(Note: This dataset was provisioned for the project but its deep integration is noted as an area for future work)*.
- **Auxiliary Annotation Files:** Supporting files for ID/Gene mapping for methylation probes (hg19) and gene annotations.

*Note: Part of raw data files were placed in `./DATA/` directory as referenced in the processing scripts.*

### Data Sources with Links
- **Gene Expression (TCGA):** [TOIL RSEM gene TPM](https://xenabrowser.net/datapages/?dataset=tcga_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net) 
- [**Gene Expression (GTEx):** [GTEx RSEM isoform TPM](https://xenabrowser.net/datapages/?dataset=gtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net) 
- **TF-Target Associations:** [ENCODE from Harmonizome](https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets) 
- **Methylation 450k (TCGA):** [PAAD](https://xenabrowser.net/datapages/?dataset-TCGA.PAAD.sampleMap%2FHuman%20Methylation450), [LUAD](https://xenabrowser.net/datapages/?dataset=TCGA.LUAD.sampleMap%2FHuman%20Methylation450), [SKCM](https://xenabrowser.net/datapages/?dataset-TCGA.SKCM.sampleMap%2FHuman%20Methylation450), [PRAD](https://xenabrowser.net/datapages/?dataset=TCGA.PRAD.sampleMap%2FHuman%20Methylation450), [STAD](https://xenabrowser.net/datapages/?dataset%20TCGA.STAD.sampleMap%2FHuman%20Methylation%20450)
- **Methylation Probe Map:** [illuminaMethyl450_hg19](https://xenabrowser.net) 
- **Gene Annotation (hg19):** Sourced via R packages `annotatr` and `TxDb.Hsapiens.UCSC.hg19.knownGene`.
- **Clinical Data (TCGA):** [Survival Supplemental Table S1](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)) 
- **Phenotypes (GTEx):** [GTEx Phenotype](https://xenabrowser.net/datapages/?dataset=GTEX_phenotype&host=https%3A%2F%2Ftoil.xenahubs.net)

## 5. Methodology
The project is executed through a sequential, three-stage analytical pipeline:

### Stage 1: Multi-Omics Data Preprocessing and Integration
This foundational stage focuses on filtering and fusing disparate data types into a cohesive, analysis-ready format.
1.  **Target Gene Identification:** The project begins by identifying all genes known to be regulated by MIR100HG from a global gene interaction database.
2.  **Expression Data Filtering:** The gene expression matrices for each cancer type are filtered to retain only MIR100HG and its identified target genes.
3.  **CpG Probe Mapping:** Genomic annotation files are used to locate the promoter region of MIR100HG. Subsequently, the DNA methylation probe map is scanned to identify all CpG probes located within this specific region.
4.  **Methylation Data Filtering:** The comprehensive methylation datasets are filtered to retain only the data from the MIR100HG-related CpG probes identified above.
5.  **Data Fusion:** The filtered expression data and filtered methylation data are integrated. The expression level of MIR100HG is appended as a new feature to the methylation matrix for each corresponding patient sample, creating the final unified dataset for downstream analysis.

### Stage 2: Correlation Analysis and Network Inference
This stage performs statistical analysis and initial network construction.
1.  **Methylation Correlation Analysis:** A Pearson correlation analysis is conducted between the methylation levels of each CpG probe in the MIR100HG promoter and the expression level of MIR100HG itself. This tests the hypothesis of epigenetic auto-regulation. Results are visualised using heatmaps.
2.  **Downstream Network Construction:** A downstream regulatory network is built by first identifying MIR100HG target genes whose expression levels are significantly correlated with MIR100HG. Then, known interactions *between* these key targets are extracted from the interaction database and used to construct and visualise a preliminary regulatory graph with `NetworkX`.

### Stage 3: GNN-based Regulatory Network Refinement
To enhance the reliability of the inferred network, this stage employs a machine learning approach.
1.  **Feature Engineering:** A global interaction graph is constructed. Node features are engineered by combining gene expression profiles with structural graph embeddings generated using the `Node2Vec` algorithm.
2.  **Ensemble GNN Model:** An ensemble of Graph Attention Network (GAT) models is implemented. The GNNs are trained on subgraphs sampled from the global network to classify and score the likelihood of true regulatory links.
3.  **Network Prediction:** The trained ensemble model is used to predict the probability of a regulatory relationship between MIR100HG and other genes.
4.  **Final Network Construction:** A final, high-confidence regulatory network is constructed by retaining only those interactions that surpass a defined probability threshold, providing a robust model of the MIR100HG regulatory landscape.

## 6. Results & Visualisation
The key outputs of the scripts include:
- **Correlation Heatmaps:** Visualising the correlation between MIR100HG promoter methylation and gene expression across different cancer types.
- **Downstream Regulatory Network Graphs:** Network diagrams illustrating the regulatory connections between MIR100HG and its most significantly correlated target genes.
- **GNN-Predicted Network:** A refined, high-confidence network graph representing the most probable regulatory interactions involving MIR100HG.

Key Findings:
- **Context-Dependent Networks:** The co-expression networks associated with MIR100HG are highly dynamic and context-specific[cite: 104]. The set of correlated genes varies significantly depending on whether MIR100HG expression is high or low, and differs dramatically between cancer and corresponding normal tissues, suggesting significant network remodeling during tumorigenesis.
- **Prognostic Significance in STAD:** A key finding emerged from the survival analysis. High expression of MIR100HG was found to be a statistically significant predictor of poorer prognosis in Stomach Adenocarcinoma (STAD) patients across all three endpoints: OS (p=0.0003), DSS (p=0.0114), and PFI. This association was not significant in the other four cancer types, highlighting a cancer-specific prognostic role.
- **Important Regulatory Features:** Feature importance analysis from the MLP model identified several key transcription factors, including **FOXP2, PBX3, TCF12, and MAFK**, as being highly predictive of MIR100HG's expression status, reinforcing its integration into core gene regulatory networks.


## 7. Installation & Usage
To replicate the analysis, follow these steps:

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/your-project-repo.git](https://github.com/your-username/your-project-repo.git)
    cd your-project-repo
    ```
2.  **Set up the environment:**
    It is highly recommended to use a virtual environment.
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```
3.  **Install dependencies:**
    A `requirements.txt` file should be created from your environment.
    ```bash
    pip install -r requirements.txt
    ```
4.  **Data Setup:**
    Download the required datasets and place them in a directory named `NTU_DATA` within the project's root folder.

5.  **Execution Workflow:**
    The scripts should be run in the following logical order:
    - **Step 1: Data Processing:** Execute `Omics Data Integration Main.py` to perform the complete data filtering and fusion pipeline.
    ```bash
    python "Omics Data Integration Main.py"
    ```
    - **Step 2: Core Analyses:** Run the analysis scripts.
    ```bash
    python regression.py  # To analyze methylation correlation
    python TF_explore.py    # To build and visualize the downstream network
    ```
    - **Step 3: GNN Modelling:** Run the graph neural network script for advanced network refinement.
    ```bash
    python RegulationNetwork.py
    ```

## 8. Conclusion & Future Work
This project successfully establishes and executes a powerful multi-omics workflow to elucidate the complex regulatory role of the lncRNA MIR100HG in cancer. The findings provide insights into both its potential auto-regulatory mechanisms via promoter methylation and its broader influence on a network of target genes. The application of a GNN model marks a significant step towards building more accurate and reliable regulatory networks from complex biological data. 
This study also confirmed that MIR100HG's regulatory function is highly dependent on its cellular context, varying by expression level, cancer type, and tissue state (normal vs. malignant). While its general involvement with key regulators is clear, its specific interactions and clinical impact are tissue- and disease-specific. The discovery of its prognostic significance in STAD provides a strong foundation for further investigation.


**Future Directions:**
- **Clinical Integration:** Integrate clinical survival data to assess the prognostic significance of the identified MIR100HG network signature.
- **Experimental Validation:** The high-confidence interactions predicted by the GNN model serve as strong candidates for experimental validation in a wet lab setting.
- **Pan-Cancer Analysis:** Extend the analysis to a broader range of cancer types to identify both common and cancer-specific regulatory patterns of MIR100HG.
- **Integration of Methylation Data:** A primary next step is the full integration of the DNA methylation data to explore how MIR100HG expression correlates with epigenetic patterns, potentially revealing key mechanisms of gene regulation.
- **Advanced Network Inference:** Employ more complex network inference algorithms (e.g., mutual information-based methods) to capture non-linear relationships and infer directionality in the regulatory network.

## 9. License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## 10. Acknowledgements 

This research was conducted as the Data Science Group Research Project for the MSc in Data Science at the University of Bristol. 

The authors wish to express their sincere gratitude to the project supervisor, **Dr. Daniel D’Andrea**, for his invaluable guidance, constructive feedback, and consistent support throughout the entire duration of this project.

Many thanks are extended to the contributing group members, **Torsa Talukdar** and **Zhanbo Wang (Mason)**, for their brillian collaboration and significant contributions to this work.

The analytical approach and scientific direction of this project were also heavily inspired by several key publications. We particularly acknowledge the foundational research by:

* **Ottaviani, S. et al. (2018).** "TGF-β induces miR-100 and miR-125b but blocks let-7a through LIN28B controlling PDAC progression." *Nature Communications*, 9, 1845. This paper provided the crucial initial link between the TGF-β pathway and MIR100HG, forming the scientific basis for this project's investigation.
* **Papoutsoglou, P. et al. (2021).** "The noncoding MIR100HG RNA enhances the autocrine function of transforming growth factor beta signaling." *Oncogene*, 40, 3748-3765. This work further elucidated the functional relationship between MIR100HG and the TGF-β signaling pathway.
* **Wu, Y., Wang, Z., et al. (2022).** "LncmiRHG-MIR100HG: A new budding star in cancer." *Frontiers in Oncology*, 12, 997532. This review provided a broad and valuable context for understanding the role of MIR100HG across various cancers.

## 11. References

1.  **ENCODE Consortium.** (2012). The ENCODE Project: Cataloging functional elements of the human genome. *Nature*, 489, 57-74.
2.  **Uszczynska-Ratajczak, B. et al.** (2018). Towards a complete map of the human long non-coding RNA transcriptome. *Nature Reviews Genetics*, 19, 535-548.
3.  **Nemeth, K. et al.** (2024). Non-coding RNAs in disease: from mechanisms to therapeutics. *Nature Reviews Genetics*, 25, 211-232.
4.  **Richardson, L. et al.** (2023). Context-dependent TGFẞ family signalling in cell fate regulation. *Nature Reviews Molecular Cell Biology*, 24, 876-894.
5.  **Ottaviani, S. et al.** (2018). TGF-β induces miR-100 and miR-125b but blocks let-7a through LIN28B controlling PDAC progression. *Nature Communications*, 9, 1845.
6.  **Jin, B. et al.** (2012). Linking DNA methyltransferases to epigenetic marks and nucleosome structure genome-wide in human tumor cells. *Cell Reports*, 2(5), 1411-1424.
7.  **Brenner, C. et al.** (2005). Myc represses transcription through recruitment of DNA methyltransferase corepressor. *EMBO Journal*, 24(2), 336-346.
8.  **Huang, W., Li, H., Yu, Q., Xiao, W. & Wang, D. O.** (2022). LncRNA-mediated DNA methylation: an emerging mechanism in cancer and beyond. *Journal of Hematology & Oncology*, 41, 100-022-02319-z.
9.  **Papoutsoglou, P. et al.** (2021). The noncoding MIR100HG RNA enhances the autocrine function of transforming growth factor beta signaling. *Oncogene*, 40, 3748-3765.
10. **Wu, Y., Wang, Z., Yu, S., Liu, D., & Sun, L.** (2022). LncmiRHG-MIR100HG: A new budding star in cancer. *Frontiers in Oncology*, 12, 997532.
11. **GTEx Consortium.** (2020). The Genotype-Tissue Expression (GTEx) Project. *Science*, 369(6509), 1318-1320.
12. **Nie, G.-H. et al.** (2017). lncRNA C22orf32-1 contributes to the tumorigenesis of nasopharyngeal carcinoma. *Oncology Letters*, 13(6), 4487-4492.
13. **Yang, M., Lu, H., Liu, J., Wu, S., Kim, P., & Zhou, X.** (2022). lncRNAfunc: a knowledgebase of lncRNA function in human cancer. *Nucleic Acids Research*, 50(D1), D1295-D1306.
14. **Harmonizome.** (2016). Harmonizome: Integrating and Mining the ENCODE Project Data. [Online]. Available: https://maayanlab.cloud/Harmonizome/.


