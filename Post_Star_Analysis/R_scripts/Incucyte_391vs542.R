library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(limma)
library(edgeR)

counts <- read.csv("./Data/pc_counts_5-20.csv", row.names = 1)
count_matrix <- as.matrix(counts)
#count_matrix <- count_matrix[rowSums(count_matrix) > 10, ]

colData <- read.csv("./Data/ColData_3_Incucyte_391vs542.csv", row.names = 1)
# 4. Subset to only NTC and FOXA2oe samples match everything and stuff
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check

keep_samples <- colData$condition %in% c("slow", "fast")
colData <- colData[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]

count_matrix <- as.matrix(counts)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

keep <- rowSums( counts(dds) >= 10 ) >= 4
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)




dds$condition <- relevel(dds$condition, ref = "slow")
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast = c("condition", "slow", "fast"))

summary(res)

#not always necessary
#res_shrink <- lfcShrink(dds, coef = "condition_FOXA2oe_vs_NTC", type = c("apeglm"))

plotMA(res)

###################################################################
resOrdered <- res[order(res$log2FoldChange), ]  

# 4. Convert to data frame
resDF <- as.data.frame(resOrdered)
resDF$Significant <- ifelse(!is.na(resDF$padj) & resDF$padj < 0.05, "Yes", "No")
# 5. Write to CSV
write.csv(resDF, file = "./incucyte_391vs542_DEgenes.csv")

resSig <- resDF[resDF$Significant == "Yes", ]
write.csv(resSig, file = "./incucyte_391vs542_DEGs_only.csv")
################################################################################

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)

# Read the CSV with gene names as row names
data <- read.csv("./incucyte_391vs542_DEGs_only.csv", row.names = 1)

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

write.csv(x = gse_df_filtered, file = "./gse_analysis_391vs542.csv")
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
ggsave("./gse_dotplot_391vs542.jpg", plot = p, width = 10, height = 7, dpi = 600)

# OR Save as PDF 
ggsave("gse_dotplot.pdf", plot = p, width = 10, height = 8)

