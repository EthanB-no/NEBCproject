library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Load organism-specific annotation package
library(org.Mm.eg.db)

# Read the CSV with gene names as row names
data <- read.csv("./5-14_results/PCA_DE_Genes.csv", row.names = 1)

# Create named logFC vector
lgfc_vector <- data$log2FoldChange
names(lgfc_vector) <- rownames(data)
lgfc_vector <- sort(lgfc_vector, decreasing = TRUE)

# Gene set enrichment analysis
gse <- gseGO(geneList = lgfc_vector,
             ont = "ALL",
             keyType = "SYMBOL",  # Change to ENSEMBL if needed
             #nPerm = 10000,
             minGSSize = 5,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "none")  # Or "none" if preferred

# Step 1: Extract GSEA result data frame
gse_df <- as.data.frame(gse)

# Step 2: Remove all rows where Description includes "bacterium"
gse_df_filtered <- subset(gse_df, !grepl("oxidoreductase", Description, ignore.case = TRUE))

write.csv(x = gse_df_filtered, file = "./5-14_results/gse_analysis_PCA_lines.csv")
# Step 3: Create a new gseaResult object
# Keep the structure of the original object but replace @result
filtered_gse <- gse
filtered_gse@result <- gse_df_filtered


# Plot example
dotplot(filtered_gse, showCategory=6, split=".sign") + facet_grid(.~.sign)

library(ggplot2)
library(enrichplot)

# Generate your plot
p <- dotplot(filtered_gse, showCategory=6, split=".sign") + facet_grid(.~.sign)

# Save as high-res PNG
ggsave("./5-14_results/gse_dotplot_PCA.jpg", plot = p, width = 10, height = 7, dpi = 600)

# OR Save as PDF (vector graphic - great for posters)
ggsave("gse_dotplot.pdf", plot = p, width = 10, height = 8)



