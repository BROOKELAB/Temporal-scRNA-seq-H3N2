# Load required packages
library('Seurat')       # For handling and analyzing single-cell RNA-seq data
library('ggplot2')      # For creating visualizations
library('dplyr')        # For data manipulation
library('data.table')   # For fast data manipulation
library('ggpubr')       # For creating publication-ready plots
library('reshape2')     # For reshaping data (e.g., wide-to-long format)

# Load single-cell RNA-seq dataset
uninfected <- readRDS("Uninfected_082924.rds")
uninfected <- subset(uninfected, subset = InfectedStatus == "NotInfected")

# Load probability data
top_prob <- read.csv("fateProbs_20240926.csv")

# Extract metadata
metadata <- uninfected@meta.data

# Subset metadata to match probability data
barcodes_to_subset <- as.vector(row.names(metadata))
top_prob <- subset(top_prob, X %in% barcodes_to_subset)

# Update metadata
metadata$cells <- rownames(metadata)
barcodes_to_subset <- as.vector(top_prob$X)
uninfected <- subset(uninfected, cells = barcodes_to_subset)
metadata <- uninfected@meta.data

# Add probability data to metadata
high_IFNL1_state <- c(top_prob$high_IFNL1_state)
metadata$high_IFNL1_state <- high_IFNL1_state
write.csv(metadata, file = "uninfected_metadata_022425.csv", row.names = TRUE)

# Create histogram of IFNL1+ probabilities
ggplot(top_prob, aes(x = high_IFNL1_state)) +
  geom_histogram(colour = "black", fill = "white") +
  scale_x_continuous(limits = c(0.2,1), name = "Probability IFNL1+", breaks = seq(0.2, 1, by = 0.2)) +
  geom_vline(aes(xintercept = mean(high_IFNL1_state, na.rm = TRUE)), linetype = "dashed", size = 1, colour = "red") +
  labs(title = "Probability IFNL1+ in all cells")

# Load gene correlation data
genes <- read.csv("20240614_corr_fateProbs_expr_v2.csv")
genes_list <- as.data.frame(genes$gene)

# Initialize results dataframe
correlation_results <- data.frame(
  Gene = character(),
  Correlation = numeric(),
  P_Value = numeric(),
  Conf_Int_Lower = numeric(),
  Conf_Int_Upper = numeric(),
  stringsAsFactors = FALSE
)

# Set up PDF output
pdf(file="correlation_plots_pearson_022425.pdf", width = 7, height = 4)
par(mfrow = c(5, 5))

# Check structure of gene list
if (!"genes$gene" %in% colnames(genes_list)) {
  stop("The data frame genes_list does not contain the column 'genes$gene'")
}

# Compute correlations
for (i in 1:nrow(genes_list)) {
  Symbol <- genes_list$`genes$gene`[i]
  if (is.null(Symbol) || Symbol == "") next
  
  if (Symbol %in% rownames(uninfected@assays$RNA$counts)) {
    expr <- uninfected@assays$RNA$counts[Symbol, ]
    metadata[[Symbol]] <- expr
    
    metadata$high_IFNL1_state <- as.numeric(as.character(metadata$high_IFNL1_state))
    metadata[[Symbol]] <- as.numeric(as.character(metadata[[Symbol]]))
    
    if (sd(metadata$high_IFNL1_state) == 0 || sd(metadata[[Symbol]]) == 0) next
    
    cor_test <- cor.test(metadata$high_IFNL1_state, metadata[[Symbol]], method = "pearson")
    correlation <- cor_test$estimate
    p_value <- cor_test$p.value
    
    conf_int_lower <- if ("conf.int" %in% names(cor_test)) cor_test$conf.int[1] else NA
    conf_int_upper <- if ("conf.int" %in% names(cor_test)) cor_test$conf.int[2] else NA
    
    if (!is.na(correlation)) {
      correlation_results <- rbind(correlation_results, data.frame(
        Gene = Symbol,
        Correlation = correlation,
        P_Value = p_value,
        Conf_Int_Lower = conf_int_lower,
        Conf_Int_Upper = conf_int_upper,
        stringsAsFactors = FALSE
      ))
      
      plot <- ggpubr::ggscatter(metadata, x = "high_IFNL1_state", y = Symbol, 
                                add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                                conf.int = TRUE, cor.coef = TRUE, cor.method = 'pearson',
                                xlab = 'Probability IFNL1+ State', ylab = 'Gene Counts', 
                                size = 3, title = paste(Symbol),
                                cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top"))
      print(plot)
    }
  }
}

dev.off()

# Save correlation results
write.csv(correlation_results, file = "correlation_results_pearson_022425.csv", row.names = FALSE)
------------------------------------------------------------------------------------------------

# Read the correlation data from a CSV file
correlation_data <- read.csv("correlation_results_pearson_022425.csv")

# Sort the data in descending order based on the correlation values
correlation_data <- correlation_data %>% arrange(desc(Correlation)) 

# Reorder the 'Gene' column factor levels according to the sorted correlation values
correlation_data$Gene <- factor(correlation_data$Gene, levels = correlation_data$Gene[order(correlation_data$Correlation,
                                                                                            decreasing= TRUE)])

# Keep only the top 20 rows (genes) with the highest correlation
correlation_data <- correlation_data[1:20,]

# Create a heatmap-like plot to visualize the correlation of the top 20 genes (Figure 3A)
ggplot(correlation_data, aes(x = Gene, y = 1, fill = Correlation)) +
  # Add a tile for each gene with a color representing the correlation value
  geom_tile(color = "black") +
  # Use a gradient color scale from white to red to represent correlation values
  scale_fill_gradient(low = "white", high = "red", limits=c(0,0.6)) +
  # Set labels for the x and y axes, as well as the title
  labs(x = "Genes", 
       y = "", 
       title = "Gene Correlations for High IFNL State (0 hrs only)") +
  # Maintain a fixed aspect ratio for better visualization
  coord_fixed(ratio = 1) +
  # Apply a classic theme for the plot
  theme_classic() +
  # Customize text and axis labels
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    # Italicize x-axis labels for better readability
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),  
    axis.text.y = element_blank(),  # Hide the y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    legend.title = element_blank(), # Remove legend title
    legend.text = element_text(size = 8, color = "black"),  # Set legend text size and color
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0.5),
    legend.key.size = unit(0.1, "in"),  # Set legend key size
    legend.position = "right",  # Place the legend to the right
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add a black border around the plot panel
    axis.line = element_blank()  # Remove axis lines
  )

# Save the plot as a PNG file with the specified dimensions and resolution
ggsave("gene_correaltion_022425.png",width = 5,height = 1.5,units = "in", dpi = 300)

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
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] reshape2_1.4.4     patchwork_1.2.0    plotly_4.10.4      ggpubr_0.6.0       data.table_1.15.4  dittoSeq_1.14.3    edgeR_4.0.16      
[8] limma_3.58.1       cowplot_1.1.3      lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
[15] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0    Seurat_5.1.0       SeuratObject_5.0.2
[22] sp_2.1-4          

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.3.2               later_1.3.2                 bitops_1.0-7                polyclip_1.10-6            
[6] fastDummies_1.7.3           lifecycle_1.0.4             rstatix_0.7.2               globals_0.16.3              lattice_0.22-6             
[11] MASS_7.3-60.0.1             backports_1.5.0             magrittr_2.0.3              rmarkdown_2.27              yaml_2.3.8                 
[16] httpuv_1.6.15               sctransform_0.4.1           spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.37.0          
[21] pbapply_1.7-2               RColorBrewer_1.1-3          abind_1.4-5                 zlibbioc_1.48.2             Rtsne_0.17                 
[26] GenomicRanges_1.54.1        BiocGenerics_0.48.1         RCurl_1.98-1.14             GenomeInfoDbData_1.2.11     IRanges_2.36.0             
[31] S4Vectors_0.40.2            ggrepel_0.9.5               irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-2       
[36] pheatmap_1.0.12             goftest_1.2-3               RSpectra_0.16-1             spatstat.random_3.2-3       fitdistrplus_1.1-11        
[41] parallelly_1.37.1           DelayedArray_0.28.0         leiden_0.4.3.1              codetools_0.2-20            tidyselect_1.2.1           
[46] matrixStats_1.3.0           stats4_4.3.2                spatstat.explore_3.2-7      jsonlite_1.8.8              progressr_0.14.0           
[51] ggridges_0.5.6              survival_3.7-0              tools_4.3.2                 ica_1.0-3                   Rcpp_1.0.12                
[56] glue_1.7.0                  SparseArray_1.2.4           gridExtra_2.3               xfun_0.45                   MatrixGenerics_1.14.0      
[61] GenomeInfoDb_1.38.8         withr_3.0.0                 BiocManager_1.30.23         fastmap_1.2.0               fansi_1.0.6                
[66] digest_0.6.36               timechange_0.3.0            R6_2.5.1                    mime_0.12                   colorspace_2.1-0           
[71] scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         utf8_1.2.4                  generics_0.1.3             
[76] httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1              uwot_0.2.2                  pkgconfig_2.0.3            
[81] gtable_0.3.5                lmtest_0.9-40               SingleCellExperiment_1.24.0 XVector_0.42.0              htmltools_0.5.8.1          
[86] carData_3.0-5               dotCall64_1.1-1             scales_1.3.0                Biobase_2.62.0              png_0.1-8                  
[91] knitr_1.47                  rstudioapi_0.16.0           tzdb_0.4.0                  nlme_3.1-165                zoo_1.8-12                 
[96] KernSmooth_2.23-24          parallel_4.3.2              miniUI_0.1.1.1              pillar_1.9.0                grid_4.3.2                 
[101] vctrs_0.6.5                 RANN_2.6.1                  promises_1.3.0              car_3.1-2                   xtable_1.8-4               
[106] cluster_2.1.6               evaluate_0.24.0             cli_3.6.3                   locfit_1.5-9.10             compiler_4.3.2             
[111] crayon_1.5.3                rlang_1.1.4                 future.apply_1.11.2         ggsignif_0.6.4              plyr_1.8.9                 
[116] stringi_1.8.4               viridisLite_0.4.2           deldir_2.0-4                munsell_0.5.1               lazyeval_0.2.2             
[121] spatstat.geom_3.2-9         Matrix_1.6-5                RcppHNSW_0.6.0              hms_1.1.3                   future_1.33.2              
[126] statmod_1.5.0               shiny_1.8.1.1               SummarizedExperiment_1.32.0 ROCR_1.0-11                 igraph_2.0.3               
[131] broom_1.0.6  

