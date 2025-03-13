library(Seurat)  # Load Seurat package for single-cell RNA-seq analysis
library(sctransform)  # Load sctransform for normalization of scRNA-seq data
library(glmGamPoi)  # Load glmGamPoi for modeling gene counts with a negative binomial distribution
library(scater)  # Load scater for single-cell RNA-seq analysis and visualization
library(scran)  # Load scran for single-cell RNA-seq analysis, especially for quality control
library(dplyr)  # Load dplyr for data manipulation and transformation
library(ggplot2)  # Load ggplot2 for creating visualizations
library(tidyr)  # Load tidyr for data tidying and reshaping
library(forcats)  # Added forcats for handling factors and working with categorical data

#BiocManager::install("hdf5r",force = TRUE)  # Install hdf5r package if needed
#BiocManager::install("simpleSingleCell")  # Install simpleSingleCell package if needed

# Read aggregate file using Read10X_h5 function
lib1 <- Read10X_h5("/outs/filtered_feature_bc_matrix.h5")
lib3 <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
lib4 <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
lib4 <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")

# Create Seurat object from the filtered feature barcode matrix for each dataset
lib1 <- CreateSeuratObject(counts = lib1, min.cells = 4, min.features = 400, project = "Ryan_scRNAseq", names.field = 2, names.delim = "-")
lib2 <- CreateSeuratObject(counts = lib2, min.cells = 4, min.features = 400, project = "Ryan_scRNAseq", names.field = 2, names.delim = "-")
lib3 <- CreateSeuratObject(counts = lib3, min.cells = 4, min.features = 400, project = "Ryan_scRNAseq", names.field = 2, names.delim = "-")
lib4 <- CreateSeuratObject(counts = lib4, min.cells = 4, min.features = 400, project = "Ryan_scRNAseq", names.field = 2, names.delim = "-")

## Subset based on mitochondria genes to filter out dead cells (remove cells with high mitochondrial content)
lib1[["percent.mt"]] <- PercentageFeatureSet(lib1, pattern = "^MT-")  # Calculate percentage of mitochondrial genes
lib1 <- subset(lib1, subset = percent.mt < 15)  # Keep cells with < 15% mitochondrial genes

lib2[["percent.mt"]] <- PercentageFeatureSet(lib2, pattern = "^MT-")
lib2 <- subset(lib2, subset = percent.mt < 25)

lib3[["percent.mt"]] <- PercentageFeatureSet(lib3, pattern = "^MT-")
lib3 <- subset(lib3, subset = percent.mt < 15)

lib4[["percent.mt"]] <- PercentageFeatureSet(lib4, pattern = "^MT-")
lib4 <- subset(lib4, subset = percent.mt < 15)

# Normalize data using the LogNormalize method
lib1 <- NormalizeData(lib1, normalization.method = "LogNormalize", scale.factor = 10000)
lib2 <- NormalizeData(lib2, normalization.method = "LogNormalize", scale.factor = 10000)
lib3 <- NormalizeData(lib3, normalization.method = "LogNormalize", scale.factor = 10000)
lib4 <- NormalizeData(lib4, normalization.method = "LogNormalize", scale.factor = 10000)

# Find the top 2000 variable genes for each dataset
lib1 <- FindVariableFeatures(lib1, selection.method = "vst", nfeatures = 2000)
lib2 <- FindVariableFeatures(lib2, selection.method = "vst", nfeatures = 2000)
lib3 <- FindVariableFeatures(lib3, selection.method = "vst", nfeatures = 2000)
lib4 <- FindVariableFeatures(lib4, selection.method = "vst", nfeatures = 2000)

# Perform cell cycle scoring using the S and G2M phase genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

lib1 <- CellCycleScoring(lib1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lib2 <- CellCycleScoring(lib2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lib3 <- CellCycleScoring(lib3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lib4 <- CellCycleScoring(lib4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Scale data after normalization, regressing out S and G2M phase scores
lib1 <- SCTransform(lib1, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = TRUE)
lib2 <- SCTransform(lib2, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = TRUE)
lib3 <- SCTransform(lib3, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = TRUE)
lib4 <- SCTransform(lib4, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = TRUE)

# Integration of multiple datasets (lib1, lib2, lib3, lib4)
lib1$library <- "Library1"
lib2$library <- "Library2"
lib3$library <- "Library3"
lib4$library <- "Library4"

# Ensure each object is a Seurat object
libraries <- list(lib1, lib2, lib3, lib4)
libraries <- libraries[sapply(libraries, inherits, "Seurat")]

# Select integration features (top 3000 variable genes)
Integration_features <- SelectIntegrationFeatures(object.list = libraries, nfeatures = 3000)

# Prepare for SCT integration
uninfected_lib <- PrepSCTIntegration(object.list = libraries, anchor.features = Integration_features)

# Find integration anchors across datasets
Integration_anchors <- FindIntegrationAnchors(object.list = uninfected_lib, normalization.method = "SCT", anchor.features = Integration_features)

# Integrate data
uninfected_lib <- IntegrateData(anchorset = Integration_anchors, normalization.method = "SCT")

# Save integrated dataset
saveRDS(uninfected_lib,"unifected_libraries_preUMAP_022625.rds")

# Run PCA
uninfected_lib <- RunPCA(uninfected_lib, verbose = FALSE)
ElbowPlot(uninfected_lib)  # Plot PCA elbow plot

# Find clusters based on PCA dimensions
uninfected_lib <- FindNeighbors(uninfected_lib, dims = 1:20)
uninfected_lib <- FindClusters(uninfected_lib, resolution = 0.1)

# Run UMAP for visualization of clusters
uninfected_lib <- RunUMAP(uninfected_lib, dims = 1:20, verbose = FALSE)

# Visualize UMAP
DimPlot(uninfected_lib, reduction = "umap")

# Join layers for RNA assay
uninfected_lib <- JoinLayers(uninfected_lib,assay = "RNA")

# Save final UMAP visualization dataset
saveRDS(uninfected_lib, "uninfected_libraries_postUMAP_022524.rds")
------------------------------------------------------------------------------------------------

# Install and load the required package for gene annotation
BiocManager::install("HGNChelper", force = TRUE)
library(HGNChelper)

# Source necessary external functions for gene set preparation and scoring
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

# Specify the database for tissue-specific gene sets
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# Define the tissue type for which to prepare gene sets (e.g., Lung)
tissue = "Lung"
gs_list = gene_sets_prepare(db_, tissue)  # Prepare the gene sets for the chosen tissue

# Load additional functions for cell-type assignment
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")

# Perform cell-type scoring using the prepared gene sets on the RNA-seq data
es.max = sctype_score(
  scRNAseqData = uninfected_lib[["integrated"]]@scale.data, 
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Assign cell types to each cluster based on scoring
cL_resutls = do.call("rbind", lapply(unique(uninfected_lib@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(uninfected_lib@meta.data[uninfected_lib@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(uninfected_lib@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# Mark low-confidence clusters as "Unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 2] = "Unknown"
sctype_scores[, 1:3]

# Assign custom classification labels to the Seurat metadata based on scoring results
uninfected_lib@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  uninfected_lib@meta.data$customclassif[uninfected_lib$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Convert custom classification to factor for easier handling
uninfected_lib@meta.data$customclassif <- as.factor(uninfected_lib@meta.data$customclassif)

# Add specific gene counts (e.g., KRT13, KRT5) to the metadata
genes_of_interest <- c("KRT13", "KRT5", "KRT14", "MKI67", "TOP2A", "TP63")
for (gene in genes_of_interest) {
  gene_counts <- uninfected@assays$RNA$count[gene, ]
  uninfected_metadata[gene] <- gene_counts
}

# Convert gene counts to long format for easier plotting
genes <- c("KRT13", "KRT5", "KRT14", "MKI67", "TOP2A", "TP63")
uninfected_long <- uninfected_metadata %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "counts")

# Define custom colors for cell types in the plot
celltype_colors <- c("#dad7cd", "#e07a5f","#3d405b","#81b29a","#f2cc8f","#9c6644")

# Create violin plot for gene expression by cell type (S6 Fig)
ggplot(uninfected_long, aes(x = customclassif, y = counts)) +
  geom_violin(aes(fill = customclassif), position = position_dodge(width = 0.9), scale = "width", bw = 0.15) + 
  scale_fill_manual(values = celltype_colors) +
  scale_y_log10() +
  labs(x = "Cell Types", y = "Gene Counts") +
  facet_wrap(~ factor(gene, levels = c("KRT13", "KRT5", "KRT14", "MKI67", "TOP2A", "TP63")), ncol = 3) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.title.x = element_blank(),  
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.position = "none",
    strip.text = element_text(face = "italic", size = 8)
  )

# Save the plot to a file
ggsave("cycling_basal.png", width = 5, height = 4, units = "in", dpi = 300)

# Update metadata with cell type names and set Seurat object identities
uninfected_lib@meta.data <- uninfected_libraries_cal@meta.data %>%
  mutate(., customclassif = fct_recode(customclassif, 
                                       "Club cells" = "Clara cells",
                                       "Secretory cells" = "Airway goblet cells",
                                       "Basal cells"= "Basal cells (Airway progenitor cells)",
                                       "Cycling basal cells" = "Immune system cells"))

# Set Seurat identities based on cell types
Idents(uninfected_lib) <- factor(Idents(uninfected_lib), 
                                 levels = c("Ciliated cells", 
                                            "Basal cells", 
                                            "Club cells", 
                                            "Ionocytes", 
                                            "Secretory cells",
                                            "Cycling basal cells"))
------------------------------------------------------------------------------------------------

# Define custom colors for cell types
celltype_colors <- c("#dad7cd", "#e07a5f","#3d405b","#81b29a","#f2cc8f","#9c6644")

# Create UMAP plot with annotated cell types (Figure 5A)
DimPlot(uninfected_lib, reduction = "umap", label = FALSE, repel = TRUE) +
  scale_color_manual(values = celltype_colors) +  # Assign colors to cell types
  labs(color = "Cell Type", title = "") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(hjust = 1),  
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "right"
  )

# Save UMAP plot to a file
ggsave("uninfecetd_annotated_UMAP_expression_libraries_022525.png", width = 4.75, height = 3.75, units = "in", dpi = 300)

# Add library metadata and recode donor information
uninfected_lib@meta.data <- uninfected_lib@meta.data %>%
  mutate(., library = fct_recode(library, "Donor 1" = "Library1", "Donor 2" = "Library2", "Donor 3"= "Library3","Donor 4" = "Library4"))

# Create UMAP plot split by library (S5 Fig)
DimPlot(uninfected_lib, reduction = "umap", label = FALSE, repel = TRUE, split.by = "library") +
  labs(color = "Cell Type", title = "") +
  scale_color_manual(values = celltype_colors) +  # Assign colors to cell types
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(hjust = 1),  
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    legend.position = "right"
  )

# Save the plot to a file
ggsave("uninfecetd_splitlibraries_022525.png", width = 4.75, height = 3.75, units = "in", dpi = 300)

# Load metadata from the Seurat object
uninfected_metadata <- uninfected_lib@meta.data

# Retrieve gene expression data for specific genes and add them to the metadata
# Add gene expression data for several genes to the uninfected metadata
genes_of_interest <- c('OASL', 'IFIT3', 'ISG15', 'IRF1', 'DDX60')
for (gene in genes_of_interest) {
  gene_counts <- uninfected@assays$RNA$count[gene, ]
  uninfected_metadata[gene] <- gene_counts
}

# List of genes for plotting and analysis
genes <- c('OASL', 'IFIT3', 'ISG15', 'IRF1', 'DDX60')

# Define custom colors for donors
library_colors <- c("Donor 1" = "#007cb5", "Donor 2" = "#d95f02",
                    "Donor 3" = "#7570b3", "Donor 4" = "#e7298a")

# Read metadata from a CSV file
uninfected_metadata <- read.csv("uninfected_metadata_new.csv")

# Pivot data for genes of interest (long format)
uninfected_long <- uninfected_metadata %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "counts"
  )

# Calculate fraction of cells expressing each gene
fraction_by_gene <- uninfected_long %>%
  group_by(library, gene) %>%
  summarise(fraction = mean(counts > 0), .groups = "drop") %>%
  complete(library, gene, fill = list(fraction = 0))  # Ensure all combinations of library and gene are present

# Recode library labels for better readability
fraction_by_gene <- fraction_by_gene %>%
  mutate(., library = fct_recode(library, "Donor 1" = "Library1", "Donor 2" = "Library2", "Donor 3" = "Library3", "Donor 4" = "Library4"))

# Create a scatter plot of the fraction of cells expressing each gene by donor (Figure 5B)
ggplot(fraction_by_gene, aes(x = factor(gene, levels = c("OASL", "IFIT3", "DDX60", "ISG15", "IRF1")), 
                             y = fraction, color = library)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +  # Fraction on y-axis
  scale_color_manual(values = library_colors) +  # Apply custom colors
  labs(x = "Genes", y = "Fraction of Cells", color = "Donor") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.title.x = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, face ="italic"),
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    strip.text = element_text(face = "italic", size = 8)  
  )

# Save the plot of fraction of cells as a PNG file
ggsave("uninfected_fraction_lib.png", width = 3.5, height = 2.5, units = "in", dpi = 300)

# Create a jitter plot of gene expression counts by donor(Figure 5C)
ggplot(uninfected_long, aes(x = library, y = counts, fill = library, color = library)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 2, dodge.width = 1), 
              alpha = 0.5, size = 0.25) +  # Adjust transparency and point size
  scale_y_log10(limits = c(1, NA)) +  # Log scale for y-axis
  scale_fill_manual(values = library_colors) +  # Set fill colors for donors
  scale_color_manual(values = library_colors) +  # Set color for legend
  labs(y = "Gene Counts", fill = "Donor", color = "Donor") +
  facet_wrap(~ factor(gene, levels = c("OASL", "IFIT3", "DDX60", "ISG15", "IRF1")), ncol = 3) +  # Facet by gene
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),  # Remove x-axis text labels
    axis.ticks.x = element_blank(),  # Remove x-axis tick marks
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),  # Make legend text bigger
    legend.key = element_blank(),  # Remove legend key box
    legend.key.size = unit(0.1, "in"),  # Increase legend dot size
    legend.position = "right",  # Position legend on the right
    strip.text = element_text(face = "italic", size = 8)  # Italicize facet labels
  )

# Save the plot as a PNG file
ggsave("uninfected_vln_expression_lib.png", width = 4.5, height = 3, units = "in", dpi = 300)

# Recode the custom class labels (cell types)
uninfected_metadata$customclassif <- factor(uninfected_metadata$customclassif, 
                                            levels = c("Ciliated cells", "Basal cells", "Club cells",
                                                       "Ionocytes", "Secretory cells", "Cycling basal cells"))

# Pivot data for genes and cell types
uninfected_long <- uninfected_metadata %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "counts"
  )

# Recode library labels again for consistency
uninfected_long <- uninfected_long %>%
  mutate(., library = fct_recode(library, "Donor 1" = "Library1", "Donor 2" = "Library2", "Donor 3" = "Library3", "Donor 4" = "Library4"))

# Count cells expressing each gene per library and cell type
cell_counts <- uninfected_long %>%
  mutate(expressed = counts > 0) %>%
  group_by(library, customclassif, gene) %>%
  summarise(count = sum(expressed), .groups = "drop")

# Compute total cells expressing each gene per library for normalization
total_cells_per_library_gene <- cell_counts %>%
  group_by(library, gene) %>%
  summarise(total = sum(count), .groups = "drop")

# Merge total counts to calculate proportions of positive cells
cell_counts <- left_join(cell_counts, total_cells_per_library_gene, by = c("library", "gene")) %>%
  mutate(fraction = count / total)

# Define custom colors for cell types
celltype_colors <- c("#dad7cd", "#e07a5f","#3d405b","#81b29a","#f2cc8f","#9c6644")

# Create stacked bar plot for the fraction of positive cells by gene and cell type (Figure 5D)
ggplot(cell_counts, aes(x = library, y = fraction, fill = customclassif)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars representing fractions
  scale_fill_manual(values = celltype_colors) +  # Use custom cell type colors
  labs(x = "Library", y = "Proportion of Positive Cells", fill = "Cell Type") +
  facet_wrap(~factor(gene, levels = c("OASL", "IFIT3", "DDX60", "ISG15", "IRF1")), ncol = 3) +  # Facet by gene
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.title.x = element_blank(),  
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),
    strip.text = element_text(face = "italic", size = 8)
  )

# Save the bar plot as a PNG file
ggsave("uninfected_bar_fract.png", width = 5, height = 3.25, units = "in", dpi = 300)

# Count the number of cells expressing each gene per library and cell type
cell_counts <- uninfected_long %>%
  mutate(expressed = counts > 0) %>%  # Create a boolean column indicating gene expression
  group_by(library, customclassif, gene) %>%  # Group by library, cell type, and gene
  summarise(count = sum(expressed), .groups = "drop")  # Sum up the number of expressing cells per group

# Compute the total number of cells per cell type per library
total_cells_per_type <- uninfected_long %>%
  group_by(library, customclassif) %>%  # Group by library and cell type
  summarise(total_cells = n(), .groups = "drop")  # Count total number of cells in each group

# Merge expression counts with total cell counts and normalize by cell type frequency
cell_counts <- left_join(cell_counts, total_cells_per_type, by = c("library", "customclassif")) %>%
  mutate(fraction_type = (count / total_cells))  # Compute the fraction of ISG-positive cells per cell type

# Define custom colors for cell types
celltype_colors <- c("#c3c0b5", "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#9c6644")

# Generate scatter plot showing the fraction of ISG-positive cells per donor (Figure 5E)
ggplot(cell_counts, aes(x = library, y = fraction_type, color = customclassif)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.6) +  # Jittered points to avoid overlap
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)) +  # Set y-axis limits and breaks
  scale_color_manual(values = celltype_colors) +  # Assign custom colors to cell types
  labs(x = "Donors", y = "Fraction ISG Positive Cells", color = "Cell Type") +  # Axis and legend labels
  theme_classic() +  # Apply a clean theme
  facet_wrap(~factor(gene, levels = c("OASL", "IFIT3", "DDX60", "ISG15", "IRF1")), ncol = 3) +  # Separate plots by gene, arranged in 3 columns
  theme(
    text = element_text(family = "Arial", color = "black"),  # Set text font and color
    axis.title = element_text(face = "bold", size = 10, color = "black"),  # Bold axis titles
    axis.title.x = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Center x-axis title
    axis.text = element_text(size = 8, color = "black"),  # Set axis text size and color
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels for readability
    legend.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),  # Bold legend title
    legend.text = element_text(size = 8, color = "black"),  # Set legend text size
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),  # Set plot title style
    legend.key.size = unit(0.1, "in"),  # Adjust legend key size
    strip.text = element_text(face = "italic", size = 8)  # Italicize facet labels
  )

# Save the plot as a PNG file
ggsave("uninfected_bar_fract_celltype.png",width = 5, height = 3.25, units = "in", dpi = 300)  # Save at 300 dpi for high resolution

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
