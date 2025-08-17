# Libraries needed for analysis
library("Seurat")       # scRNA-seq analysis framework
library('tidyverse')    # data manipulation and visualization (ggplot2, dplyr, etc.)
library('ggplot2')      # plotting
library("cowplot")      # combining plots
library("edgeR")        # differential expression tools
library('dplyr')        # data wrangling
library('data.table')   # fast data frame manipulation
library('ggpubr')       # publication-ready plots
library('devtools')     # package development tools
library('msigdbr')        # MSigDB gene sets (hallmarks, C2, C5, etc.)
library('dplyr')          # loaded again (redundant but harmless)



# Load Seurat object and probability data
allcells_postfilt <- readRDS("data_seurat_postUMAP_081824.rds")   # main Seurat object
prob <- read.csv("fateProbs_20240926.csv")                        # cell fate probabilities

# Set the correct column as rownames so they match Seurat cells
rownames(prob) <- prob$cells   # assumes "cells" column contains cell barcodes
prob$cells <- NULL             # remove redundant column

# Align Seurat object with probabilities
common_cells <- intersect(colnames(allcells_postfilt), rownames(prob))   # overlapping cells

# Subset Seurat object and prob table to common cells
allcells_postfilt <- subset(allcells_postfilt, cells = common_cells)
fate_probs <- prob[common_cells, , drop = FALSE]

# Add fate probabilities to Seurat metadata
allcells_postfilt <- AddMetaData(allcells_postfilt, metadata = fate_probs)

# Extract metadata
meta <- allcells_postfilt@meta.data

# Set up probability thresholds
threshold <- 0.3
prob_columns <- c("high_IFNL1_state", "low_IFNL1_state_1", 
                  "low_IFNL1_state_2", "low_IFNL1_state_3")
probs <- meta[, prob_columns]   # select fate probability columns

# Find markers
allcells_postfilt <- readRDS("allcells_postfilt_endstate.rds")   # reload object with end state annotation

# Differential expression between "end_state" groups
markers <- FindAllMarkers(
  allcells_postfilt,
  group.by = "end_state",   # cluster/grouping variable
  only.pos = FALSE          # include both up- and down-regulated genes
)

# Save results
saveRDS(markers, file = "marker_072025.rds")
saveRDS(allcells_postfilt, "allcells_postfilt_endstate.rds")

# Load hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")

gene_categories <- hallmark_sets %>%
  dplyr::select(gene_symbol, gs_name) %>%
  dplyr::distinct() %>%
  dplyr::rename(Gene = gene_symbol, Category = gs_name)

# Expression matrix setup
expr_mat <- GetAssayData(allcells_postfilt, slot = "scale.data")   # use scaled expression

# Ensure metadata column is character
allcells_postfilt$end_state <- as.character(allcells_postfilt$end_state)
meta <- allcells_postfilt@meta.data

# Convert expression to long format for joining with metadata
library('tidyr')
library('tibble')

expr_long <- expr_mat %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
  left_join(meta[, "end_state", drop = FALSE] %>% rownames_to_column("Cell"), by = "Cell")

# Join with gene sets & summarise
expr_cat <- expr_long %>%
  inner_join(gene_categories, by = "Gene") %>%
  group_by(Category, end_state) %>%
  summarise(AvgExpr = mean(Expression, na.rm = TRUE), .groups = "drop")

# Reshape into matrix for heatmap
heatmap_matrix <- expr_cat %>%
  pivot_wider(names_from = end_state, values_from = AvgExpr) %>%
  column_to_rownames("Category") %>%
  as.matrix()

# Filter gene categories by variability
heatmap_matrix <- readRDS("heatmap_matrix.rds")   # reload saved matrix

library(matrixStats)

row_vars <- rowVars(as.matrix(heatmap_matrix))   # row-wise variance
row_means <- rowMeans(heatmap_matrix)
cv <- row_vars / row_means                       # coefficient of variation
top_cv <- names(sort(cv, decreasing = TRUE))[1:25]   # select top 25 variable categories
heatmap_matrix_filtered <- heatmap_matrix[top_cv, ]

# Clean up row/col names
rownames(heatmap_matrix_filtered) <- rownames(heatmap_matrix_filtered) %>%
  gsub("^HALLMARK_", "", .) %>%   # remove "HALLMARK_" prefix
  gsub("_", " ", .)               # replace underscores with spaces

colnames(heatmap_matrix_filtered) <- c("High IFNL", "Low IFNL 1", "Low IFNL 2", "Low IFNL 3")

# Plot heatmap (all states) S3 Fig A
library('ComplexHeatmap')
library('circlize')
library('svglite')

# Scale rows (z-score per gene set)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_filtered)))

ht <- Heatmap(
  heatmap_matrix_scaled,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 8, fontfamily = "Arial")),
  col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "darkred")),
  column_title = "High/Low IFNL States",
  column_title_gp = gpar(fontfamily = "Arial", fontsize = 12, just = "right"),
  row_names_side = "right",
  show_row_dend = FALSE
)

# Save as SVG
svglite("/Users/joelrivera/Downloads/heatmap_output_all.svg", width = 4, height = 6, system_fonts = list(sans = "Arial"))
draw(ht)
dev.off()



# Repeat: Heatmap for only Low IFNL states  S3 Fig B
heatmap_matrix <- readRDS("/Users/joelrivera/Downloads/heatmap_matrix.rds")
heatmap_matrix <- heatmap_matrix[, !(colnames(heatmap_matrix) %in% "high_IFNL1")]   # drop high IFNL

# Recalculate top variable gene sets
row_vars <- rowVars(as.matrix(heatmap_matrix))
row_means <- rowMeans(heatmap_matrix)
cv <- row_vars / row_means
top_cv <- names(sort(cv, decreasing = TRUE))[1:25]
heatmap_matrix_filtered <- heatmap_matrix[top_cv, ]

rownames(heatmap_matrix_filtered) <- rownames(heatmap_matrix_filtered) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("_", " ", .)

colnames(heatmap_matrix_filtered) <- c("Low IFNL 1", "Low IFNL 2", "Low IFNL 3")

# Scale & plot again
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_filtered)))

ht <- Heatmap(
  heatmap_matrix_scaled,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 8, fontfamily = "Arial")),
  col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "darkred")),
  column_title = "Low IFNL States",
  column_title_gp = gpar(fontfamily = "Arial", fontsize = 12, just = "right"),
  row_names_side = "right",
  show_row_dend = FALSE
)

svglite("/Users/joelrivera/Downloads/heatmap_output_low.svg", width = 4, height = 6, system_fonts = list(sans = "Arial"))
draw(ht)
dev.off()

─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────
setting  value
version  R version 4.3.2 (2023-10-31)
os       macOS 15.6
system   x86_64, darwin20
ui       RStudio
language (EN)
collate  en_US.UTF-8
ctype    en_US.UTF-8
tz       America/New_York
date     2025-08-17
rstudio  2024.04.1+748 Chocolate Cosmos (desktop)
pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/ (via rmarkdown)
quarto   1.4.553 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────
package          * version    date (UTC) lib source
abind              1.4-8      2024-09-12 [1] CRAN (R 4.3.3)
assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.3.0)
babelgene          22.9       2022-09-29 [1] CRAN (R 4.3.0)
backports          1.5.0      2024-05-23 [1] CRAN (R 4.3.3)
BiocGenerics       0.48.1     2023-11-01 [1] Bioconductor
broom              1.0.8      2025-03-28 [1] CRAN (R 4.3.3)
cachem             1.1.0      2024-05-16 [1] CRAN (R 4.3.3)
car                3.1-3      2024-09-27 [1] CRAN (R 4.3.3)
carData            3.0-5      2022-01-06 [1] CRAN (R 4.3.0)
circlize         * 0.4.16     2024-02-20 [1] CRAN (R 4.3.2)
cli                3.6.5      2025-04-23 [1] CRAN (R 4.3.3)
clue               0.3-66     2024-11-13 [1] CRAN (R 4.3.3)
cluster            2.1.8.1    2025-03-12 [1] CRAN (R 4.3.3)
codetools          0.2-20     2024-03-31 [1] CRAN (R 4.3.2)
colorspace         2.1-1      2024-07-26 [1] CRAN (R 4.3.3)
ComplexHeatmap   * 2.18.0     2023-10-24 [1] Bioconductor
cowplot          * 1.2.0      2025-07-07 [1] CRAN (R 4.3.3)
crayon             1.5.3      2024-06-20 [1] CRAN (R 4.3.3)
curl               6.4.0      2025-06-22 [1] CRAN (R 4.3.3)
data.table       * 1.17.8     2025-07-10 [1] CRAN (R 4.3.3)
deldir             2.0-4      2024-02-28 [1] CRAN (R 4.3.2)
devtools         * 2.4.5      2022-10-11 [1] CRAN (R 4.3.0)
dichromat          2.0-0.1    2022-05-02 [1] CRAN (R 4.3.0)
digest             0.6.37     2024-08-19 [1] CRAN (R 4.3.3)
doParallel         1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
dotCall64          1.2        2024-10-04 [1] CRAN (R 4.3.3)
dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.3.0)
edgeR            * 4.0.16     2024-02-18 [1] Bioconductor 3.18 (R 4.3.2)
ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
evaluate           1.0.4      2025-06-18 [1] CRAN (R 4.3.3)
farver             2.1.2      2024-05-13 [1] CRAN (R 4.3.3)
fastDummies        1.7.5      2025-01-20 [1] CRAN (R 4.3.3)
fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.3.3)
fitdistrplus       1.2-4      2025-07-03 [1] CRAN (R 4.3.3)
forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
foreach            1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
Formula            1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
fs                 1.6.6      2025-04-12 [1] CRAN (R 4.3.3)
future             1.58.0     2025-06-05 [1] CRAN (R 4.3.2)
future.apply       1.20.0     2025-06-06 [1] CRAN (R 4.3.2)
generics           0.1.4      2025-05-09 [1] CRAN (R 4.3.3)
GetoptLong         1.0.5      2020-12-15 [1] CRAN (R 4.3.0)
ggplot2          * 3.5.2      2025-04-09 [1] CRAN (R 4.3.3)
ggpubr           * 0.6.1      2025-06-27 [1] CRAN (R 4.3.3)
ggrepel            0.9.6      2024-09-07 [1] CRAN (R 4.3.3)
ggridges           0.5.6      2024-01-23 [1] CRAN (R 4.3.2)
ggsignif           0.6.4      2022-10-13 [1] CRAN (R 4.3.0)
GlobalOptions      0.1.2      2020-06-10 [1] CRAN (R 4.3.0)
globals            0.18.0     2025-05-08 [1] CRAN (R 4.3.2)
glue               1.8.0      2024-09-30 [1] CRAN (R 4.3.3)
goftest            1.2-3      2021-10-07 [1] CRAN (R 4.3.0)
gridExtra          2.3        2017-09-09 [1] CRAN (R 4.3.0)
gtable             0.3.6      2024-10-25 [1] CRAN (R 4.3.3)
hms                1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.3.2)
htmlwidgets        1.6.4      2023-12-06 [1] CRAN (R 4.3.0)
httpuv             1.6.16     2025-04-16 [1] CRAN (R 4.3.3)
httr               1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
ica                1.0-3      2022-07-08 [1] CRAN (R 4.3.0)
igraph             2.1.4      2025-01-23 [1] CRAN (R 4.3.3)
IRanges            2.36.0     2023-10-24 [1] Bioconductor
irlba              2.3.5.1    2022-10-03 [1] CRAN (R 4.3.0)
iterators          1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
jsonlite           2.0.0      2025-03-27 [1] CRAN (R 4.3.3)
KernSmooth         2.23-26    2025-01-01 [1] CRAN (R 4.3.3)
knitr              1.50       2025-03-16 [1] CRAN (R 4.3.3)
later              1.4.2      2025-04-08 [1] CRAN (R 4.3.3)
lattice            0.22-7     2025-04-02 [1] CRAN (R 4.3.3)
lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.3.0)
limma            * 3.58.1     2023-10-31 [1] Bioconductor
listenv            0.9.1      2024-01-29 [1] CRAN (R 4.3.2)
lmtest             0.9-40     2022-03-21 [1] CRAN (R 4.3.0)
locfit             1.5-9.12   2025-03-05 [1] CRAN (R 4.3.3)
lubridate        * 1.9.4      2024-12-08 [1] CRAN (R 4.3.3)
magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
MASS               7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.0)
Matrix             1.6-5      2024-01-11 [1] CRAN (R 4.3.0)
matrixStats        1.5.0      2025-01-07 [1] CRAN (R 4.3.3)
memoise            2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
mime               0.13       2025-03-17 [1] CRAN (R 4.3.3)
miniUI             0.1.2      2025-04-17 [1] CRAN (R 4.3.3)
msigdbr          * 25.1.1     2025-07-21 [1] CRAN (R 4.3.3)
nlme               3.1-168    2025-03-31 [1] CRAN (R 4.3.3)
parallelly         1.45.0     2025-06-02 [1] CRAN (R 4.3.3)
patchwork          1.3.1      2025-06-21 [1] CRAN (R 4.3.3)
pbapply            1.7-4      2025-07-20 [1] CRAN (R 4.3.3)
pillar             1.11.0     2025-07-04 [1] CRAN (R 4.3.3)
pkgbuild           1.4.8      2025-05-26 [1] CRAN (R 4.3.3)
pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
pkgload            1.4.0      2024-06-28 [1] CRAN (R 4.3.3)
plotly             4.11.0     2025-06-19 [1] CRAN (R 4.3.3)
plyr               1.8.9      2023-10-02 [1] CRAN (R 4.3.0)
png                0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
polyclip           1.10-7     2024-07-23 [1] CRAN (R 4.3.3)
profvis            0.4.0      2024-09-20 [1] CRAN (R 4.3.3)
progressr          0.15.1     2024-11-22 [1] CRAN (R 4.3.3)
promises           1.3.3      2025-05-29 [1] CRAN (R 4.3.3)
purrr            * 1.1.0      2025-07-10 [1] CRAN (R 4.3.3)
R6                 2.6.1      2025-02-15 [1] CRAN (R 4.3.3)
RANN               2.6.2      2024-08-25 [1] CRAN (R 4.3.3)
RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
Rcpp               1.1.0      2025-07-02 [1] CRAN (R 4.3.3)
RcppAnnoy          0.0.22     2024-01-23 [1] CRAN (R 4.3.2)
RcppHNSW           0.6.0      2024-02-04 [1] CRAN (R 4.3.2)
readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.3.0)
remotes            2.5.0      2024-03-17 [1] CRAN (R 4.3.2)
reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
reticulate         1.43.0     2025-07-21 [1] CRAN (R 4.3.3)
rjson              0.2.23     2024-09-16 [1] CRAN (R 4.3.3)
rlang              1.1.6      2025-04-11 [1] CRAN (R 4.3.3)
rmarkdown          2.29       2024-11-04 [1] CRAN (R 4.3.3)
ROCR               1.0-11     2020-05-02 [1] CRAN (R 4.3.0)
RSpectra           0.16-2     2024-07-18 [1] CRAN (R 4.3.3)
rstatix            0.7.2      2023-02-01 [1] CRAN (R 4.3.0)
rstudioapi         0.17.1     2024-10-22 [1] CRAN (R 4.3.3)
Rtsne              0.17       2023-12-07 [1] CRAN (R 4.3.0)
S4Vectors          0.40.2     2023-11-23 [1] Bioconductor
scales             1.4.0      2025-04-24 [1] CRAN (R 4.3.3)
scattermore        1.2        2023-06-12 [1] CRAN (R 4.3.0)
sctransform        0.4.2      2025-04-30 [1] CRAN (R 4.3.3)
sessioninfo        1.2.3      2025-02-05 [1] CRAN (R 4.3.3)
Seurat           * 5.3.0      2025-04-23 [1] CRAN (R 4.3.3)
SeuratObject     * 5.1.0      2025-04-22 [1] CRAN (R 4.3.3)
shape              1.4.6.1    2024-02-23 [1] CRAN (R 4.3.2)
shiny              1.11.1     2025-07-03 [1] CRAN (R 4.3.3)
sp               * 2.2-0      2025-02-01 [1] CRAN (R 4.3.3)
spam               2.11-1     2025-01-20 [1] CRAN (R 4.3.3)
spatstat.data      3.1-6      2025-03-17 [1] CRAN (R 4.3.3)
spatstat.explore   3.5-2      2025-07-22 [1] CRAN (R 4.3.3)
spatstat.geom      3.5-0      2025-07-20 [1] CRAN (R 4.3.3)
spatstat.random    3.4-1      2025-05-20 [1] CRAN (R 4.3.3)
spatstat.sparse    3.1-0      2024-06-21 [1] CRAN (R 4.3.3)
spatstat.univar    3.1-4      2025-07-13 [1] CRAN (R 4.3.3)
spatstat.utils     3.1-5      2025-07-17 [1] CRAN (R 4.3.3)
statmod            1.5.0      2023-01-06 [1] CRAN (R 4.3.0)
stringi            1.8.7      2025-03-27 [1] CRAN (R 4.3.3)
stringr          * 1.5.1      2023-11-14 [1] CRAN (R 4.3.0)
survival           3.8-3      2024-12-17 [1] CRAN (R 4.3.3)
svglite          * 2.2.1      2025-05-12 [1] CRAN (R 4.3.3)
systemfonts        1.2.3      2025-04-30 [1] CRAN (R 4.3.3)
tensor             1.5.1      2025-06-17 [1] CRAN (R 4.3.3)
textshaping        1.0.1      2025-05-01 [1] CRAN (R 4.3.3)
tibble           * 3.3.0      2025-06-08 [1] CRAN (R 4.3.3)
tidyr            * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.3.2)
tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
timechange         0.3.0      2024-01-18 [1] CRAN (R 4.3.0)
tzdb               0.5.0      2025-03-15 [1] CRAN (R 4.3.3)
urlchecker         1.0.1      2021-11-30 [1] CRAN (R 4.3.0)
usethis          * 3.1.0      2024-11-26 [1] CRAN (R 4.3.3)
uwot               0.2.3      2025-02-24 [1] CRAN (R 4.3.3)
vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.3.0)
viridisLite        0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
withr              3.0.2      2024-10-28 [1] CRAN (R 4.3.3)
xfun               0.52       2025-04-02 [1] CRAN (R 4.3.3)
xtable             1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
yaml               2.3.10     2024-07-26 [1] CRAN (R 4.3.3)
zoo                1.8-14     2025-04-10 [1] CRAN (R 4.3.3)
