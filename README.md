# **Temporal-scRNA-seq-H3N2**
### **10X Chromium Single-Cell RNA-seq Analysis**  
#### **Intrinsic OASL expression licenses interferon induction during influenza A virus infection**

## **Overview**  
This repository contains scripts and analysis pipelines used to quantify the expression of interferons and interferon-stimulated genes (ISGs) in infected and uninfected cells using single-cell RNA sequencing (scRNA-seq).  

The reference assembly scripts can be found in **BROOKELAB/SingleCell**. The analysis was performed using **Cell Ranger** on a **10X Chromium Single-Cell** dataset, including:  

- Human alveolar epithelial cells (**A549**) infected with **A/Perth/16/2009 (Perth09)** at different time points.  
- Uninfected/untreated human bronchial epithelial cells (**HBECs**).  

### **Analysis Pipeline Includes**  
✅ Quality control and filtering  
✅ Normalization and annotation  
✅ Dimensionality reduction  
✅ Quantification of gene expression frequencies  
✅ Correlation analysis  

---

## **Requirements**  

Ensure that **R (v4.3.2)** is installed and accessible from the command line.  

### **Required R Packages**  
The following R packages are required:  
- `Seurat (v5.1.0)`  
- `scater (v1.30.0)`  
- `scran (v1.30.2)`  
- `sctransform (v0.4.1)`  
- `ggplot2 (v3.5.1)`  
- `dplyr (v1.1.4)`  
- `glmGamPoi (v1.12.2)`  
- `tidyr (v1.3.1)`  
- `forcats (v1.0.0)`  

---

## **Preliminary Processing for A549 Cells**  

Raw reads were demultiplexed and mapped to a host-virus hybrid reference using the **10X Chromium Cell Ranger** software:  
🔗 [10X Genomics Single-Cell Gene Expression Software](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)  

### **Steps**  

1️⃣ **Cell Ranger Count**  
- Aligns reads, filters cells, and performs barcode/UMI counting.  
- Input: FASTQ files  
- Output: CSV, BAM, MEX, H5 files  
- Input file: `files_A549.txt`  
- Script: `cellranger_count_A549.sh`  

2️⃣ **Cell Ranger Aggregation (aggr)**  
- Aggregates samples across experimental conditions and replicates.  
- Input: `.h5` files  
- Output: TSV files and filtered matrices for downstream analysis  
- Input file: `AggrList_A549.csv`  
- Script: `aggr_all_A549.sh`  

---

## **A549 scRNA-seq Analysis**  

### **Script: `A549_seurat.R`**  

This script performs:  
✅ Preliminary filtering (empty drops, cell cycle calling, doublet removal, gene/cell filtering)  
✅ Infection status determination based on viral read percentages  
✅ Normalization and creation of a **Seurat object** for downstream analysis  

**Gene Correlation Analysis:**  
- Uses outputs from `Gene_Correlation.R` to compute correlation coefficients between gene expression and the probability of transitioning to a **high IFNL state**.  
- Input:  
  - Filtered **Seurat** object (mock/0-hour timepoint).  
  - Data frames with transition probabilities per cell.  
  - Gene list for correlation analysis.  
- Output:  
  - **CSV file**: Correlation coefficients, p-values, and confidence intervals.  
  - **PDF file**: Correlation plots.  

Additional analyses:  
✅ Temporal-NoSpliceVelo transition probability analysis  
✅ Filtering of genes based on correlation coefficients  
✅ Differential expression analysis in mock populations  

---

## **Preliminary Processing for HBECs**  

Raw reads were demultiplexed and mapped to a **host-virus hybrid reference** using **10X Chromium Cell Ranger**.  
For validation, sequences were aligned against the **H1N1 A/California/07/2009 genome** to confirm the mock populations.  

🔗 [10X Genomics Single-Cell Gene Expression Software](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)  

### **Steps**  

1️⃣ **Cell Ranger Count**  
- Aligns reads, filters cells, and performs barcode/UMI counting.  
- Input: FASTQ files  
- Output: CSV, BAM, MEX, H5 files  
- Input file: `files_HBEC.txt`  
- Script: `cellranger_count_HBEC.sh`  

---

## **HBEC scRNA-seq Analysis**  

### **Script: `HBEC_seurat.R`**  

This script performs:  
✅ Preliminary filtering (empty drops, cell cycle calling, doublet removal, gene/cell filtering)  
✅ Integration of different scRNA-seq libraries into one **Seurat object**
✅ Normalization of a **Seurat object** for downstream analysis  

Additional analyses:  
✅ **Cell type annotation** using [ScType](https://github.com/IanevskiAleksandr/sc-type)  
✅ **Manual annotation** for misclassified cells based on gene expression profiles  
✅ **Intrinsic ISG expression assessment** across different cell types  
