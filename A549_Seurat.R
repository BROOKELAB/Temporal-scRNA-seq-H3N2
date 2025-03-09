# Libraries needed for analysis
library(dittoSeq)         # For advanced plotting for Seurat objects
library(tidyverse)        # General data manipulation and visualization
library(Seurat)           # For single-cell RNA sequencing analysis
library(sctransform)      # For normalization of scRNA-seq data
library(glmGamPoi)        # For glmGamPoi method for variance stabilization in scRNA-seq
library(simpleSingleCell) # For preprocessing of single-cell RNA-seq data
library(scater)           # For scRNA-seq quality control and analysis
library(scran)            # For analysis of single-cell RNA-seq data
library(magrittr)         # For piping (%>%)
library(edgeR)            # For differential expression analysis of RNA-seq data
library(dplyr)            # For data manipulation
library(cowplot)          # For advanced plotting functionality

# Loading the data (from HDF5 file)
all_data <- Read10X_h5("filtered_feature_bc_matrix.h5")

# Create Seurat object and process the data
allcells_prefilt <- CreateSeuratObject(counts = all_data, min.cells = 4,
                                       min.features = 400,
                                       project = "Perth09_timecourse", names.field = 2, names.delim = "-")

# Read sample information and process it
sample_order <- read.csv("outs/aggregation.csv")
sample_order$sample_id[1:11] <- c("Uninfected", "Perth09_8hrs_1", "Perth09_8hrs_2", "Perth09_10hrs_1", "Perth09_10hrs_2", 
                                  "Perth09_12hrs_1", "Perth09_12hrs_2", "Perth09_14hrs_1", "Perth09_14hrs_2", 
                                  "Perth09_16hrs_1", "Perth09_16hrs_2")
sample_order$number <- as.character(1:nrow(sample_order))
levels(allcells_prefilt$orig.ident) <- sample_order$sample_id[order(sample_order$number)]
allcells_prefilt$orig.ident <- factor(allcells_prefilt$orig.ident, levels = sample_order$sample_id)

# Add metadata from aggregation csv
rownames(sample_order) <- sample_order$sample_id
allcells_prefilt$time <- sample_order[allcells_prefilt$orig.ident, "time"]

# Quality control for mitochondrial genes (remove dead cells)
allcells_prefilt[["percent.mt"]] <- PercentageFeatureSet(allcells_prefilt, pattern = "^MT-")
allcells_prefilt <- subset(allcells_prefilt, subset = percent.mt < 15)

# QC visualization
VlnPlot(allcells_prefilt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("data_mtQC_081724.pdf")

# Further QC visualization with FeatureScatter
plot1 <- FeatureScatter(allcells_prefilt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(allcells_prefilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("data_featureQC_081724.pdf", plot = plot1 + plot2)

# Save data after QC
saveRDS(allcells_prefilt, "data_seurat_preCellCycle_081724.rds")
allcells_prefilt <- readRDS("data_seurat_preCellCycle_081724.rds")

# Add viral gene counts and calculate viral percentage
segment_gene_names <- c("NS", "NA.", "HA", "NP", "PB2", "PB1", "PA", "M")
temp <- GetAssayData(allcells_prefilt, slot = "count")
allcells_prefilt$host_counts <- colSums(temp[!rownames(temp) %in% segment_gene_names,])
allcells_prefilt$virus_counts <- colSums(temp[segment_gene_names,])
allcells_prefilt$virus_pct <- allcells_prefilt$virus_counts / allcells_prefilt$nCount_RNA * 100

# Log transformation of virus percentage
halfmin_virus <- min(allcells_prefilt$virus_pct[allcells_prefilt$virus_pct > 0]) / 2
allcells_prefilt$virus_pct_log10 <- log10(allcells_prefilt$virus_pct + halfmin_virus)

# Identification of infected cells based on viral gene expression
des.all <- density(allcells_prefilt$virus_pct_log10[allcells_prefilt$virus_counts > 0])
min.all <- des.all$x[which(diff(sign(diff(des.all$y))) == 2) + 1]
min.all.use <- -0.8
allcells_prefilt$InfectedStatus <- ifelse(allcells_prefilt$virus_pct_log10 > min.all.use, "Infected", "NotInfected")

# Add viral gene metadata to Seurat object
temp <- GetAssayData(allcells_prefilt, slot = "count")[segment_gene_names,] %>% as.matrix() %>% t() %>% as.data.frame()
allcells_prefilt <- AddMetaData(allcells_prefilt, temp)
allcells_prefilt <- allcells_prefilt[!rownames(allcells_prefilt) %in% segment_gene_names]

# Normalize data and find variable genes
allcells_prefilt <- NormalizeData(allcells_prefilt, normalization.method = "LogNormalize", scale.factor = 10000)
allcells_prefilt <- FindVariableFeatures(allcells_prefilt, selection.method = "vst", nfeatures = 2000)

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
allcells_prefilt <- CellCycleScoring(allcells_prefilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Scale the data and run PCA
allcells_prefilt <- SCTransform(allcells_prefilt, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = TRUE)
allcells_prefilt <- RunPCA(allcells_prefilt, verbose = FALSE)

# Find clusters and UMAP visualization
allcells_prefilt <- FindNeighbors(allcells_prefilt, dims = 1:20)
allcells_prefilt <- FindClusters(allcells_prefilt, resolution = 0.3)
allcells_prefilt <- RunUMAP(allcells_prefilt, dims = 1:40, verbose = FALSE)

# Save the post-UMAP Seurat object
saveRDS(allcells_prefilt, "data_seurat_postUMAP_081824.rds")
------------------------------------------------------------------------------------------------

# Recode time variable for better readability
allcells_postfilt@meta.data <- allcells_postfilt@meta.data %>%
  mutate(., time = fct_recode(time, "0 hrs" = "0hrs", "8 hrs" = "8hrs", "10 hrs" = "10hrs", "12 hrs" = "12hrs",
                              "14 hrs" = "14hrs", "16 hrs" = "16hrs"))

# Set the 'time' variable as the identity class for cells
Idents(allcells_postfilt) <- "time"

# Ensure the identity levels are ordered by the time points
Idents(allcells_postfilt) <- factor(Idents(allcells_postfilt), 
                                    levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs"))

# Define color palette for the time points
time_colors <- c("#84847F", "#ecb500", "limegreen", "darkturquoise", "steelblue2", "orchid3")

# Plot UMAP, colored by time points, using the defined colors (Figure 1B)
DimPlot(allcells_postfilt, reduction = "umap", label = FALSE, repel = TRUE, alpha = 0.3) +
  scale_color_manual(values = time_colors) +  # Assign custom colors to time points
  labs(color = "Time", title = "") +
  theme_classic() +  # Apply a classic theme to the plot
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(hjust = 1),  # Adjust x-axis label alignment
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "right"
  )

ggsave(filename = 'split_time_030225.png', width = 4.75, height = 4.25, units = "in", dpi = 300)


# Identify cells with IFNL1 expression greater than 0
infl_pos <- WhichCells(allcells_postfilt, expression = IFNL1 > 0)

# Recode the 'time' variable to have a specific factor level order
allcells_postfilt$time <- factor(allcells_postfilt$time, levels = c("0hrs", "8hrs", "10hrs", "12hrs", "14hrs", "16hrs"))

# DimPlot to visualize UMAP with the highlighted cells having IFNL1 expression (Figure 1C)
DimPlot(allcells_postfilt, 
        split.by = "time",  # Split plots by time points
        cells.highlight = infl_pos,  # Highlight cells with IFNL1 expression
        cols.highlight = "red",  # Color for highlighted cells
        cols = "grey",  # Color for non-highlighted cells
        pt.size = 0.2,  # Adjust size of non-highlighted points
        sizes.highlight = 0.2,  # Adjust size of highlighted points
        ncol = 3)  # Number of columns for the layout of plots

ggsave(filename = 'IFNL1_infection_status_011425.jpeg', height = 15, width = 20, units = "cm",dpi = 600)

# Extract and add gene expression counts for several genes to the metadata
metadata <- allcells_postfilt@meta.data
genes_of_interest <- c('IFNL1', 'IFNB1','IFNL2','MX1', 'ISG15', 'TRIM25')
for (gene in genes_of_interest) {
  gene_counts <- allcells_postfilt@assays$RNA$count[gene, ]
  metadata[gene] <- gene_counts
}

# Save the updated metadata to a CSV file
write.csv(metadata, "metadata_all_021825.csv")

# Read the saved metadata from CSV
metadata <- read_csv('metadata_all_021825.csv')

# Recode time variable for the imported metadata
metadata <- metadata %>%
  mutate(., time = fct_recode(time, "0 hrs" = "0hrs", "8 hrs" = "8hrs", "10 hrs" = "10hrs", "12 hrs" = "12hrs",
                              "14 hrs" = "14hrs", "16 hrs" = "16hrs"))

# Calculate the fraction of cells expressing IFNL1 across time points
fraction_by_time <- metadata %>%
  group_by(time) %>%
  summarise(fraction = mean(IFNL1 > 0))

# Create a bar plot showing the fraction of cells expressing IFNL1 (Figure 1D)
ggplot(fraction_by_time, aes(x = factor(time, levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs")),
                             y = fraction)) +
  geom_bar(stat = "identity", fill = "#84847F", color = "black") +  # Bar color and borders
  scale_y_continuous(limits = c(0, 0.06)) +  # Set y-axis limits
  labs(x = "Times", 
       y = expression(atop(bold("Fraction of ")*bolditalic(IFNL1) * bold(" Cells")))) +  # Label formatting
  theme_classic() +  # Apply a classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "none"
  )

ggsave("fraction_IFNL1_021825.png", width = 2.6, height = 2, units = "in", dpi = 300)

# Create a violin plot of IFNL1 gene expression across time points (Figure 1E)
ggplot(metadata, aes(x = factor(time, levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs")),
                     y = IFNL1)) +
  geom_violin(fill = "#84847F", bw = 0.2, adjust = 0.8) +  # Create violin plot
  geom_jitter(width = 0.2, color = "black", size = .5) +  # Add jittered points for each cell
  scale_y_log10() +  # Log scale for better visualization of expression
  labs(x = "Times", 
       y = expression(bolditalic(IFNL1) * bold(" Counts"))) +  # Axis labels
  theme_classic() +  # Classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "none"
  )

ggsave("counts_IFNL1_021825.png", width = 2.6, height = 2, units = "in", dpi = 300)
------------------------------------------------------------------------------------------------

# Read metadata CSV file containing cell information
metadata <- read_csv('metadata_all_021825.csv')

# Filter out samples where time is 0hrs
metadata_ifn <- metadata[metadata$time != "0hrs", ]

# Define the list of genes to analyze
genes <- c("IFNB1", "IFNL2")  # Replace with your genes of interest

# Compute the fraction of expressing cells for each gene across different time points
fraction_by_time <- map_dfr(genes, function(gene) {
  metadata_ifn %>%
    group_by(time) %>%
    summarise(fraction = mean(.data[[gene]] > 0), .groups = "drop") %>%
    mutate(Gene = gene)  # Add a column for the gene name
})

# Plot the fraction of expressing cells for each gene across time points (S1 Fig A)
ggplot(fraction_by_time, aes(x = factor(time, levels = c("8hrs", "10hrs", "12hrs", "14hrs", "16hrs")),
                             y = fraction, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +  # Create bar plot
  scale_y_continuous(limits = c(0, 0.06)) +  # Limit the y-axis from 0 to 0.06
  labs(x = "Time", y = "Fraction of Expressing Cells", fill = "Gene") +  # Add labels
  scale_fill_manual(values = c("#84847F","#264653")) +  # Customize color palette
  theme_classic() +  # Apply classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),  # Set font
    axis.title = element_text(face = "bold", size = 10, color = "black"),  # Bold axis titles
    axis.text = element_text(size = 8, color = "black"),  # Set axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),  # Bold legend title
    legend.text = element_text(size = 8, color = "black",face = "italic"),  # Italicize legend text
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Center plot title
    legend.key.size = unit(0.1, "in")  # Set legend key size
  )

# Save the fraction plot as a PNG image
ggsave("fraction_IFNs_022125.png", width = 3, height = 2, units = "in", dpi = 300)

# Reshape data to long format for the violin plot
metadata_long <- metadata_ifn %>%
  select(time, all_of(genes)) %>%  # Select time and gene expression columns
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")  # Reshape data to long format

# Plot the gene expression distribution (violins) for each gene over time (S1 Fig C)
ggplot(metadata_long, aes(x = factor(time, levels = c("0hrs", "8hrs", "10hrs", "12hrs", "14hrs", "16hrs")),
                          y = Expression, fill = Gene)) +
  geom_violin(position = position_dodge(width = 0.9), scale = "width", trim = TRUE, bw = 0.2) +  # Create violin plot
  geom_jitter(aes(color = Gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), 
              size = 0.25, alpha = 0.5) +  # Add jittered points to the plot
  scale_y_log10() +  # Log-transform the y-axis
  scale_fill_manual(values = c("#84847F","#264653")) +  # Customize color palette
  scale_color_manual(values = c("black","black"), guide = "none") +  # Set color for jitter points
  labs(x = "Times", y = "Gene Counts", fill = "Gene") +  # Add labels
  theme_classic() +  # Apply classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),  # Set font
    axis.title = element_text(face = "bold", size = 10, color = "black"),  # Bold axis titles
    axis.text = element_text(size = 8, color = "black"),  # Set axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),  # Bold legend title
    legend.text = element_text(size = 8, color = "black", face = "italic"),  # Italicize legend text
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Center plot title
    legend.key.size = unit(0.1, "in")  # Set legend key size
  )

# Save the violin plot as a PNG image
ggsave("counts_IFNs_022125.png", width = 3, height = 2, units = "in", dpi = 300)

# Define another list of genes for ISGs (MX1, TRIM25, ISG15)
genes <- c("MX1","TRIM25", "ISG15")  # Replace with your genes of interest

# Compute the fraction of expressing cells for the new set of genes
fraction_by_time <- map_dfr(genes, function(gene) {
  metadata_ifn %>%
    group_by(time) %>%
    summarise(fraction = mean(.data[[gene]] > 0), .groups = "drop") %>%
    mutate(Gene = gene)  # Add a column for the gene name
})

# Plot the fraction of expressing cells for the new set of genes across time points (S1 Fig B)
ggplot(fraction_by_time, aes(x = factor(time, levels = c("8hrs", "10hrs", "12hrs", "14hrs", "16hrs")),
                             y = fraction, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +  # Create bar plot
  scale_y_continuous(limits = c(0, 0.5)) +  # Limit the y-axis from 0 to 0.5
  labs(x = "Time", y = "Fraction of Expressing Cells", fill = "Gene") +  # Add labels
  scale_fill_manual(values = c("#84847F","#264653", "#f4a261")) +  # Customize color palette
  theme_classic() +  # Apply classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),  # Set font
    axis.title = element_text(face = "bold", size = 10, color = "black"),  # Bold axis titles
    axis.text = element_text(size = 8, color = "black"),  # Set axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),  # Bold legend title
    legend.text = element_text(size = 8, color = "black", face = "italic"),  # Italicize legend text
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Center plot title
    legend.key.size = unit(0.1, "in")  # Set legend key size
  )

# Save the fraction plot as a PNG image
ggsave("fraction_ISGs_022125.png", width = 3.4, height = 2, units = "in", dpi = 300)

# Reshape data to long format for the violin plot
metadata_long <- metadata_ifn %>%
  select(time, all_of(genes)) %>%  # Select time and gene expression columns
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")  # Reshape data to long format

# Plot the gene expression distribution (violins) for the ISGs over time (S1 Fig D)
ggplot(metadata_long, aes(x = factor(time, levels = c("8hrs", "10hrs", "12hrs", "14hrs", "16hrs")),
                          y = Expression, fill = Gene)) +
  geom_violin(position = position_dodge(width = 0.9), scale = "width", trim = TRUE, bw = 0.2) +  # Create violin plot
  geom_jitter(aes(color = Gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), 
              size = 0.25, alpha = 0.5) +  # Add jittered points to the plot
  scale_y_log10() +  # Log-transform the y-axis
  scale_fill_manual(values = c("#84847F","#264653", "#f4a261")) +  # Customize color palette
  scale_color_manual(values = c("black","black", "black"), guide = "none") +  # Set color for jitter points
  labs(x = "Times", y = "Gene Counts", fill = "Gene") +  # Add labels
  theme_classic() +  # Apply classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),  # Set font
    axis.title = element_text(face = "bold", size = 10, color = "black"),  # Bold axis titles
    axis.text = element_text(size = 8, color = "black"),  # Set axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),  # Bold legend title
    legend.text = element_text(size = 8, color = "black", face = "italic"),  # Italicize legend text
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Center plot title
    legend.key.size = unit(0.1, "in")  # Set legend key size
  )

# Save the violin plot as a PNG image
ggsave("counts_ISGs_022125.png", width = 3.4, height = 2, units = "in", dpi = 300)
------------------------------------------------------------------------------------------------

# Read the CSV file containing immune probabilities for cells
all_immune_prob <- read.csv("fateProbs_20240926.csv")

# Merge the immune probabilities data with metadata based on a common column ("cells")
merged_df <- full_join(all_immune_prob, metadata, by="cells")

# Remove rows with missing values (NA)
merged_df <- na.omit(merged_df)

# Plot the density distribution of high IFNL state probabilities for each time point (Figure 2C)
ggplot(merged_df, aes(x = high_IFNL1_state, 
                      colour = factor(time, levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs")))) +
  # Create a density plot with transparency (alpha), specific bandwidth (bw), and line width
  geom_density(alpha = 0.2, bw = 0.02, linewidth = 0.75) +  
  theme_minimal() +  # Use minimal theme for a clean look
  scale_y_log10(limits = c(0.000001, 1000000)) +  # Log scale for y-axis with custom limits
  scale_x_log10() +  # Log scale for x-axis
  # Custom color scale for the time points
  scale_color_manual(values = c("#84847F", "#ecb500", "limegreen", "darkturquoise", "steelblue2", "orchid3")) +  
  # Title and labels for the axes
  labs(x = "Transition Probabilities for \nHigh IFNL State", y = "Density", title = "", colour = "Time") +  
  theme_classic() +  # Classic theme with white background
  # Customize text and axis labels
  theme(
    text = element_text(family = "Arial", color = "black"),               # Global text color and font
    axis.title = element_text(face = "bold", size = 10, color = "black"),    # Axis titles
    axis.text = element_text(size = 8, color = "black"),                   # Axis tick labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),  # Legend title
    legend.text = element_text(size = 8, color = "black"),                 # Legend text
    plot.title = element_text(face = "bold", size = 10, color = "black",hjust = 0.5),    # Plot title
    legend.key.size = unit(0.1, "in"),                                     # Reduce the size of legend keys
    legend.position = "right"  # Place the legend to the right
  )

# Save the plot as a PNG file with specific dimensions and resolution
ggsave("probability_density_high_022425.png",width = 3,height = 2,units = "in", dpi = 300)

# Calculate the fraction of cells with high IFNL1 state probability > 0.3 for each time point
fraction_by_time <- merged_df %>%
  group_by(time) %>%
  summarise(fraction = mean(high_IFNL1_state > 0.3))  # Fraction of cells with high IFNL1 state

# Plot the fraction of cells with high IFNL1 state probability > 0.3 by time point (Figure 2D)
ggplot(fraction_by_time, aes(x = factor(time, levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs")),
                             y = fraction,
                             fill = factor(time, levels = c("0 hrs", "8 hrs", "10 hrs", "12 hrs", "14 hrs", "16 hrs")))) +
  # Create a bar plot for the fraction of cells per time point
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 0.1)) +  # Set y-axis limits
  scale_fill_manual(values = c("#84847F", "#ecb500", "limegreen", "darkturquoise", "steelblue2", "orchid3")) +  # Custom color scale
  labs(x = "Times", y = "Fraction of Cells with \nProbability >0.3",
       title = "", fill = "Times") +  # Labels and title
  theme_classic() +  # Apply classic theme
  # Customize the appearance of the plot
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black",hjust = 0.5),
    legend.key.size = unit(0.1, "in"),  # Reduce legend key size
    legend.position = "none"  # Remove the legend
  )

# Save the fraction plot as a PNG file with specific dimensions and resolution
ggsave("fraction_high_022425.png",width = 2.5,height = 2,units = "in", dpi = 300)
------------------------------------------------------------------------------------------------

# Subset the uninfected cells and extract metadata
uninfected <- subset(allcells_postfilt, subset = time == "0hrs")
uninfected_metadata <- uninfected@meta.data

# Add gene expression data for several genes to the uninfected metadata
genes_of_interest <- c('OASL', 'IFIT3', 'ISG15', 'IRF1', 'DDX60')
for (gene in genes_of_interest) {
  gene_counts <- uninfected@assays$RNA$count[gene, ]
  uninfected_metadata[gene] <- gene_counts
}

# Save the uninfected metadata to CSV
write.csv(uninfected_metadata, "metadata_uninfected_021825.csv")

# Create a long-format dataset for the uninfected cells and gene counts
uninfected_metadata_long <- uninfected_metadata %>%
  pivot_longer(
    cols = all_of(genes_of_interest),
    names_to = "gene",
    values_to = "counts"
  )

# Calculate the fraction of cells expressing each gene for the uninfected cells
fraction_by_gene <- uninfected_metadata_long %>%
  group_by(gene) %>%
  summarise(fraction = mean(counts > 0), .groups = "drop")

# Create a bar plot of the fraction of cells expressing each gene in uninfected cells (Figure 3C)
ggplot(fraction_by_gene, aes(x = factor(gene, levels = genes_of_interest), 
                             y = fraction)) +
  geom_bar(stat = "identity", fill = "#84847F") +  # Bar color
  scale_y_continuous(limits = c(0, 1)) +  # Limit y-axis from 0 to 1
  labs(x = "Genes", y = "Fraction of Cells") +  # Axis labels
  theme_classic() +  # Classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "none"
  )

# Save the bar plot of gene fractions (commented out)
ggsave("fraction_genes_021825.png", width = 2.6, height = 2, units = "in", dpi = 300)

# Create a violin plot of gene expression counts in uninfected cells (Figure 3D)
ggplot(uninfected_long, aes(x = factor(gene, levels = genes_of_interest), y = counts)) +
  geom_violin(fill = "#84847F") +  # Create violin plot
  geom_jitter(width = 0.2, color = "black", size = 0.5) +  # Add jittered points for each cell
  scale_y_log10() +  # Log scale for better visualization of expression
  labs(x = "Genes", y = "Gene Counts") +  # Axis labels
  theme_classic() +  # Classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "none"
  )

# Save the violin plot of gene counts (commented out)
ggsave("counts_genes_021825.png", width = 2.6, height = 2, units = "in", dpi = 300)


sessionInfo()

R version 4.3.2 (2023-10-31)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS 15.3.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] magrittr_2.0.3              scran_1.30.2                scater_1.30.1               scuttle_1.12.0              SingleCellExperiment_1.24.0
[6] SummarizedExperiment_1.32.0 Biobase_2.62.0              GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
[11] S4Vectors_0.40.2            BiocGenerics_0.48.1         MatrixGenerics_1.14.0       matrixStats_1.3.0           simpleSingleCell_1.26.0    
[16] glmGamPoi_1.12.2            sctransform_0.4.1           reshape2_1.4.4              patchwork_1.2.0             plotly_4.10.4              
[21] ggpubr_0.6.0                data.table_1.15.4           dittoSeq_1.14.3             edgeR_4.0.16                limma_3.58.1               
[26] cowplot_1.1.3               lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
[31] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1              
[36] tidyverse_2.0.0             Seurat_5.1.0                SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.3.2             later_1.3.2               bitops_1.0-7              CodeDepends_0.6.6        
[6] polyclip_1.10-6           graph_1.80.0              XML_3.99-0.17             fastDummies_1.7.3         lifecycle_1.0.4          
[11] rstatix_0.7.2             processx_3.8.4            globals_0.16.3            lattice_0.22-6            MASS_7.3-60.0.1          
[16] backports_1.5.0           rmarkdown_2.27            yaml_2.3.8                metapod_1.10.1            httpuv_1.6.15            
[21] spam_2.10-0               spatstat.sparse_3.1-0     reticulate_1.37.0         pbapply_1.7-2             RColorBrewer_1.1-3       
[26] abind_1.4-5               zlibbioc_1.48.2           Rtsne_0.17                RCurl_1.98-1.14           GenomeInfoDbData_1.2.11  
[31] ggrepel_0.9.5             irlba_2.3.5.1             listenv_0.9.1             spatstat.utils_3.1-2      pheatmap_1.0.12          
[36] goftest_1.2-3             RSpectra_0.16-1           dqrng_0.4.1               spatstat.random_3.2-3     fitdistrplus_1.1-11      
[41] parallelly_1.37.1         DelayedMatrixStats_1.24.0 DelayedArray_0.28.0       leiden_0.4.3.1            codetools_0.2-20         
[46] tidyselect_1.2.1          viridis_0.6.5             ScaledMatrix_1.10.0       spatstat.explore_3.2-7    jsonlite_1.8.8           
[51] BiocNeighbors_1.20.2      progressr_0.14.0          ggridges_0.5.6            survival_3.7-0            tools_4.3.2              
[56] ica_1.0-3                 Rcpp_1.0.12               glue_1.7.0                SparseArray_1.2.4         gridExtra_2.3            
[61] xfun_0.51                 withr_3.0.0               BiocManager_1.30.23       fastmap_1.2.0             bluster_1.12.0           
[66] fansi_1.0.6               rsvd_1.0.5                callr_3.7.6               digest_0.6.36             timechange_0.3.0         
[71] R6_2.5.1                  mime_0.12                 colorspace_2.1-0          scattermore_1.2           tensor_1.5               
[76] spatstat.data_3.1-2       utf8_1.2.4                generics_0.1.3            httr_1.4.7                htmlwidgets_1.6.4        
[81] S4Arrays_1.2.1            uwot_0.2.2                pkgconfig_2.0.3           gtable_0.3.5              lmtest_0.9-40            
[86] XVector_0.42.0            htmltools_0.5.8.1         carData_3.0-5             dotCall64_1.1-1           scales_1.3.0             
[91] png_0.1-8                 knitr_1.47                rstudioapi_0.16.0         tzdb_0.4.0                nlme_3.1-165             
[96] zoo_1.8-12                KernSmooth_2.23-24        vipor_0.4.7               parallel_4.3.2            miniUI_0.1.1.1           
[101] pillar_1.9.0              grid_4.3.2                vctrs_0.6.5               RANN_2.6.1                promises_1.3.0           
[106] BiocSingular_1.18.0       car_3.1-2                 beachmat_2.18.1           xtable_1.8-4              cluster_2.1.6            
[111] beeswarm_0.4.0            evaluate_0.24.0           cli_3.6.3                 locfit_1.5-9.10           compiler_4.3.2           
[116] crayon_1.5.3              rlang_1.1.4               future.apply_1.11.2       ggsignif_0.6.4            ps_1.7.7                 
[121] ggbeeswarm_0.7.2          plyr_1.8.9                stringi_1.8.4             BiocParallel_1.36.0       viridisLite_0.4.2        
[126] deldir_2.0-4              munsell_0.5.1             lazyeval_0.2.2            spatstat.geom_3.2-9       Matrix_1.6-5             
[131] RcppHNSW_0.6.0            hms_1.1.3                 sparseMatrixStats_1.14.0  future_1.33.2             statmod_1.5.0            
[136] shiny_1.8.1.1             ROCR_1.0-11               igraph_2.0.3              broom_1.0.6          