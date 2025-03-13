# Temporal-scRNA-seq-H3N2
### 10X Chromium single cell scripts
#### Intrinsic OASL expression licenses interferon induction during influenza A virus infection

## Analysis Overview:
This analysis are describe the approach used to quantify expression of interferona and interferon 
stimulated genes in infected and uninfected cells using scRNA-seq. Scripts required for reference 
assembly can be found in BROOKELAB/SingleCell. 

This analysis was performed using Cell Ranger on a 10X Chromium Single Cell experiment of human 
alveolar epithelial cells (A549) collected at different times during infection with A/Perth/16/2009
(Perth09). Also to Single Cell experiments of uninfected/untreated human bronchial epithelial 
cells (HBECs). The analysis performed here include quality control, filtering, normalization, annotation, 
dimensional reduction, quantification of gene expression frequencies, and correlation analysis.

### Requirements:
The scripts requires R to be installed and made available from the command line * R (v4.3.2)

The following R packages are also required: * Seurat (v5.1.0) * scater (v1.30.0) * 
scran (v1.30.2) * sctransform (v0.4.1) * ggplot2 (v3.5.1) * dplyr (v1.1.4) * 
glmGamPoi (v1.12.2) * tidyr (v1.3.1) * forcats (v1.0.0)


### Preliminary steps for A549s: 
Raw reads were demultiplexed and mapped to a host/virus hybrid reference using the 10X Chromium
Cell Ranger software package:
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

-1- Use cellranger count (cellranger v7.1.0) for alignment, filtering, and barcode/UMI counting. 
This takes uses FASTQ as input and provides multiple file outputs (CSV, BAM, MEX, and H5). 
	- txt file for input: files_A549.txt or 
	- Script for cellranger counts: cellranger_count_A549.sh

-2- Perform cellranger aggr to aggregate all samples generated from one experiment, including
experimental conditions and replicates. This analysis requires .h5 files as input and provides a
list of tsv files and filtered matrixes that can be used for downstream analysis. 
	- csv file for input: AggrList_A549.csv 
	- Script for cellranger aggr: aggr_all_A549.sh


### Script: A549_seurat.R

The script performs preliminary filtering (empty drops, cell cycle calling, filter cells by
detected genes, filter genes by min cells, double calling) and then calculates viral read
percentages to call infection status and viral gene presence/absence, resulting in normalized seurat
object for future analysis. This analysis is used to determine expression of interferons and interferon
stimulated genes (ISGs) at different times during influenza infection.

In addition this analysis uses the outputs of the script Gene_Correlation.R which calculates correlation
coefficients between cells probability to transition to high_IFNL_state an the expression of close
to 2000 genes for the mock (0 hrs) timepoint.
	- Input for the analysis is a filtered seurat object containing only mock (0 hrs) timepoint. In 
	addtion it required two data frames containing the transitionprobability calculated for each cell 
	present in the seurat object and a the list of genes to compare. The transition probabilities
	were calulcated using temporal-NoSpliceVelo.
	- Output is a csv and pdf file containing the correaltion coefficient, p-value, and confidence 
	intervals for each gene in the list and the corresponding correlation plots respectively.

This script further analyzes the transition probabilities obtained from temporal-NoSpliceVelo and
determines the probability across cells at different times during infection as well as the fraction
cells with high probabilities of transitioning into terminal state of interest. In addition it uses the
output from the Gene_Correlation.R to filter genes basede on correlation coefficient and determines 
their distinct expression within the mock (0 hrs) population.


### Preliminary steps for HBECs: 
Raw reads were demultiplexed and mapped to a host/virus hybrid reference using the 10X Chromium.
These were aligned against a reference sequence which included the H1N1 A/California/07/2009 genome 
to confirm the mock populations since some of the ones used originated from experiments comparing 
H1N1-infected and uninfected populations. 
Cell Ranger software package:
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

Use cellranger count (cellranger v7.1.0) for alignment, filtering, and barcode/UMI counting. 
This takes uses FASTQ as input and provides multiple file outputs (CSV, BAM, MEX, and H5). 
	- txt file for input: files_HBEC.txt or 
	- Script for cellranger counts: cellranger_count_HBEC.sh


### Script: HBEC_seurat.R

The script performs preliminary filtering (empty drops, cell cycle calling, filter cells by
detected genes, filter genes by min cells, double calling) and then calculates viral read
percentages to call infection status and viral gene presence/absence, resulting in normalized seurat
object for future analysis. This analysis is used to determine expression of interferons and interferon
stimulated genes (ISGs) in publicly available datasets from uninfected/untreated primary human 
bronchial epithelial cells.

In addition this script determines cell types present in the different libraries of the primary HBECs. ScType
(https://github.com/IanevskiAleksandr/sc-type) was used to determine the distinct cell types present in the 
populations. The annotation was validated and cells types that did not fit they're description were annotated 
based on expression profiles. Once cells were annotated intrinsic ISG expression was assessed across libraries
and cell types. 


