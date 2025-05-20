library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(limma)
library(edgeR)

counts <- read.csv("./5-14_results/documentation/PCA_Lines_Counts.csv", row.names = 1)
count_matrix <- as.matrix(counts)
#count_matrix <- count_matrix[rowSums(count_matrix) > 10, ]

colData <- read.csv("./5-14_results/documentation/PCA_BC_lines_Coldata.CSV", row.names = 1)
# 4. Subset to only NTC and FOXA2oe samples match everything and stuff
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check

keep_samples <- colData$condition %in% c("Group1", "Group2")
colData <- colData[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]

count_matrix <- as.matrix(counts)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

#keep <- rowSums( counts(dds) >= 10 ) >= 4
#dds <- dds[keep,]
dds <- estimateSizeFactors(dds)

#normalized_counts <- counts(dds, normalized = TRUE)
#write.csv(as.data.frame(normalized_counts), file = "normalized_counts.csv")


dds$condition <- relevel(dds$condition, ref = "Group2")
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast = c("condition", "Group1", "Group2"))

summary(res)


res_shrink <- lfcShrink(dds, coef = "condition_FOXA2oe_vs_NTC", type = c("apeglm"))

plotMA(res)

###################################################################
resOrdered <- res[order(res$padj), ]  

# 4. Convert to data frame
resDF <- as.data.frame(resOrdered)
resDF$Significant <- ifelse(!is.na(resDF$padj) & resDF$padj < 0.05, "Yes", "No")
# 5. Write to CSV
write.csv(resDF, file = "./PCA_groups2_DE_res.csv")

################################################################################

library(Glimma)

library(DESeq2)
library(Glimma)
library(edgeR)

# 1. Extract counts and convert to DGEList
counts <- counts(dds, normalized = TRUE)
dge <- DGEList(counts = counts)
dge$samples <- as.data.frame(colData(dds))  # Sample metadata

# 2. Optional: annotate genes
anno <- data.frame(GeneID = rownames(res))

# 3. Run glimmaVolcano
glimmaVolcano(
  dge = dge,
  res = res,
  status = res$padj < 0.05,
  anno = anno,
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  path = "glimma_custom",
  file = "foxa1_volcano"
)


glimmaVolcano(dds)

################################################################################
library(ggplot2)
library(ggrepel)

# Assuming your DESeq2 results object is named `res`
res$gene <- rownames(res)

# Add DE status based on padj and log2FC
res$DE_Status <- "Not DE"
res$DE_Status[res$padj < 0.05 & res$log2FoldChange > 1] <- "Over"
res$DE_Status[res$padj < 0.05 & res$log2FoldChange < -1] <- "Under"
res$DE_Status <- factor(res$DE_Status, levels = c("Under", "Not DE", "Over"))

# Colors for volcano plot
glimma_colors <- c("Under" = "#5B8DB8", "Not DE" = "grey70", "Over" = "#D33F49")

# Genes to label
genes_of_interest <- list(
  downreg <- c(
    "Foxa1", "Hoxd3", "Gria3", "Tnn", "Mfap2", "Sorbs2", "Cpa6", "Sox2",
    "Gm41741", "Ptn", "Sfrp1", "Dner", "Pcdh17", "Asic3", "Lepr", "Hoxd4",
    "Gm38462", "Hoxc10", "Meg3", "Hoxb9", "Tnfrsf11b", "Slc2a10", "Wnt10a", "Il18", "Foxp1
", "Foxg1", "Foxs1"
  ),
  upreg_2 <- c(
    "Hsd11b2", "S100g", "Gata4", "Slc16a12", "Cracdl", "E330013P04Rik", "Tmem266",
    "Hpgd", "F13a1", "Ccdc3", "Greb1", "Trpc6", "P2ry14", "Hottip", "Spon1",
    "Tnnc1", "Tnfaip8l3", "Bdkrb2", "Car2", "Ptprq", "Oasl1", "Ptgdr2", "Dach1", "Dmrta1"
  )
  
)

# Plot
jpeg("volcano_plot_PCA.jpg", width = 10, height = 10, units = "in", res = 500)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = DE_Status), alpha = 0.6) +
  scale_color_manual(values = glimma_colors) +
  geom_text_repel(
    data = subset(res, gene %in% unlist(genes_of_interest)),
    aes(label = gene),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    force = 1,
    force_pull = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(legend.position = "right")
dev.off()

dds <- DESeq(dds, quiet=TRUE)

glimmaMDS(dds)
glimmaVolcano(dds)
glimmaMA(dds)

library(Glimma)

# Convert DESeqDataSet to counts matrix and extract sample groups
counts <- counts(dds, normalized = TRUE)
group <- dds$condition

# Extract DESeq2 results
res <- results(dds)

# Run glimmaMA with all arguments
glimmaVolcano(dds, groups = group)
