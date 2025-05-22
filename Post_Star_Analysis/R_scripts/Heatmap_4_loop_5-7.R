# Load required libraries
library(DESeq2)
library(pheatmap)
library(SummarizedExperiment)
library(dplyr)
library(RColorBrewer)

# Read in data
counts <- read.csv("./Data/pc_counts_5-20.csv", row.names = 1)
count_matrix <- as.matrix(counts)

colData <- read.csv("./WDL_ColData.CSV", row.names = 1)
colData <- colData[match(colnames(count_matrix), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(count_matrix)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

# (Optional) filter low-expressed genes
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
vst <- varianceStabilizingTransformation(dds)
vst_mat <- assay(vst)

# Define sample subset (optional)
samples_to_include <- c("BC01", "BC02", "BC03", "BC04","BC05","BC07", "BC08", "BC09", "BC10")
vst_mat <- vst_mat[, samples_to_include]

vst_mat_df <- as.data.frame(vst_mat)

write.csv(vst_mat_df, "vst_biomarkers_5-16.csv")

dists <- dist(t(assay(vst)))
plot(hclust(dists))

# Define marker genes by group
marker_lists <- list(
  Luminal = c("Krt20", "Pparg", "Foxa1", "Gata3", "Snx31", "Upk1a", "Upk2", "Fgfr3"),
  ECM_SmoothMuscle = c("Pgm5", "Des", "C7", "Sfrp4", "Comp", "Sgcd"),
  EMT_Claudin = c("Zeb1", "Zeb2", "Snai1", "Twist1", "Cdh2", "Cldn3", "Cldn4", "Cldn7"),
  Basal = c("Cd44", "Krt6a", "Krt5", "Krt14", "Col17a1"),
  Squamous = c("Dsc3", "Gsdmc", "Tgm1", "Pi3", "Tp63"),
  Immune = c("Cd274", "Pdc1lg2", "Ido1", "Cxcl11", "L1cam", "Saa1"),
  Neuronal = c("Msi1", "Plekhg4b", "Gnag", "Peg10", "Rnd2", "Alpl", "Sox2", "Tubb2b", "Foxq1", "Fgfr3", "Syp", "Chga", "Cd56", "Eno2", "Ncam1" ),
  CIS = c("Crtac1", "Ctse", "Padi3", "Msn", "Nr3c1"),
  SHH = c("Shh", "Bmp5")
)

# Create output directory
output_dir <- "/home/ethan/Desktop/BioInformatics/5-14_results/Heatmaps5-14"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Color palette
brighter_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate and save one heatmap per marker group
for (group in names(marker_lists)) {
  genes <- marker_lists[[group]]
  genes_present <- genes[genes %in% rownames(vst_mat)]
  
  if (length(genes_present) == 0) {
    message("Skipping ", group, ": no genes found.")
    next
  }
  
  # Subset VST matrix to relevant genes
  vst_subset <- vst_mat[genes_present, , drop = FALSE]
  
  # Remove constant rows to avoid NA in scaling
  vst_subset <- vst_subset[apply(vst_subset, 1, sd) > 0, , drop = FALSE]
  
  # Scale gene expression 
  vst_scaled_subset <- t(scale(t(vst_subset)))
  
  # Optional: keep input gene order
  row_order <- intersect(genes, rownames(vst_scaled_subset))
  vst_scaled_subset <- vst_scaled_subset[row_order, , drop = FALSE]
  
  # Save heatmap
  jpeg_filename <- file.path(output_dir, paste0("heatmap_", group, ".jpeg"))
  jpeg(jpeg_filename, width = 2200, height = 2000, res = 400)
  pheatmap(vst_scaled_subset,
           color = brighter_colors,
           breaks = seq(-1.5, 1.5, length.out = 101),
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize_row = 16,
           fontsize_col = 10,
           show_rownames = TRUE,
           border_color = NA,
           treeheight_row = 30,
           treeheight_col = 30,
           main = paste("Marker Expression -", group))
  dev.off()
}


