# 1. Load libraries
library(limma)
library(Glimma)
library(edgeR)
library(EnhancedVolcano)
library(dplyr)

# 2. Read counts and sample info
counts <- read.csv("./Data/pc_counts_5-20.csv", row.names = 1)
colData <- read.csv("./5-14_results/documentation/PCA_BC_lines_Coldata.CSV", row.names = 1)


# 3. Reorder colData to match count matrix columns
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check


# 4. Subset to only NTC and FOXA2oe samples
keep_samples <- colData$condition %in% c("Group1", "Group2")
colData <- colData[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]


# 5. Create DGEList and filter low-expressed genes
#counts <- counts[rowSums(counts) > 10, ]
#counts are counts
dge <- DGEList(counts = counts)

#set NTC as reference for comparisons
group <- relevel(factor(colData$condition), ref = "Group2")
dge$samples$group <- group
design <- model.matrix(~ group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = TRUE]


#cpm_data <- cpm(dge)
#keep <- rowSums(cpm_data > 2) >= (ncol(dge) / 2)  # must be >1 CPM in ≥ half the samples
#dge <- dge[keep, , keep.lib.sizes = TRUE]


dge <- calcNormFactors(dge, method = "TMM")


dge$samples$lib.size            # raw total count
dge$samples$norm.factors        # scaling factor


# Run after you’ve defined dge
plotMDS(dge, labels = colData$condition, col = as.numeric(factor(colData$condition)))

# 6. Voom transformation
# 6. Voom transformation with sample quality weights
v <- voom(dge, design, plot = TRUE)

# 7. Fit linear model and empirical Bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = 2, number = Inf, sort.by = "logFC")
colnames(fit$coefficients)
res$symbol <- rownames(res)
res$Significant <- ifelse(res$adj.P.Val < 0.01 & abs(res$logFC) > 2, "Yes", "No")



# Define genes you want to label
# change these as needed
res$gene <- rownames(res)
# Add differential expression status
res$DE_Status <- "Not DE"
res$DE_Status[res$adj.P.Val < 0.05 & res$logFC > 1] <- "Over"
res$DE_Status[res$adj.P.Val < 0.05 & res$logFC < -1] <- "Under"
res$DE_Status <- factor(res$DE_Status, levels = c("Under", "Not DE", "Over"))

write.csv(res, file = "./5-13_PCA_limma.csv", row.names = TRUE)

glimma_colors <- c("Under" = "#5B8DB8", "Not DE" = "grey70", "Over" = "#D33F49")

# Define specific genes to label
genes_of_interest <- list(
  Top_Upregulated = c(
    "P2ry4", "Gm41620", "Gm36185", "Gvin-ps2", 
    "Miat", "Cx3cr1", "Hspa1b", "Acta1", "Foxa1"
  ),
  
  Top_Downregulated = c(
    "Gjc3", "Zfp213", "Ascl1", "Depp1", "Fcrlb",
    "Leng9", "Homez", "Heyl", "Tril", "Mmd2"
  )
)

jpeg("volcano_plot_glimma_style_foxa1_513.jpg", width = 9, height = 5, units = "in", res = 500)
ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
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
#title = "Differentially Expressed Genes",
x = "Log2 Fold Change",
y = "-Log10 Adjusted P-value"
) +
theme(legend.position = "right")
#plot.title = element_text(hjust = 0.5)
dev.off()

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Load organism-specific annotation package
library(org.Mm.eg.db)

# Read the CSV with gene names as row names
data <- read.csv("./5-13_results/5-13_bladder_cell.csv", row.names = 1)

# Create named logFC vector
lgfc_vector <- data$logFC
names(lgfc_vector) <- rownames(data)
lgfc_vector <- sort(lgfc_vector, decreasing = TRUE)

# Gene set enrichment analysis
gse <- gseGO(geneList = lgfc_vector,
             ont = "ALL",
             keyType = "SYMBOL",  # Change to ENSEMBL if needed
             #nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "none")  # Or "none" if preferred

# Step 1: Extract GSEA result data frame
gse_df <- as.data.frame(gse)

# Step 2: Remove all rows where Description includes "bacterium"
gse_df_filtered <- subset(gse_df, !grepl("bacterium", Description, ignore.case = TRUE))

write.csv(x = gse_df_filtered, file = "gse_analysis_bc_lines.csv")
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
ggsave("gse_dotplot_bladdercell.jpeg", plot = p, width = 10, height = 5, dpi = 600)

# OR Save as PDF (vector graphic - great for posters)
ggsave("gse_dotplot.pdf", plot = p, width = 10, height = 8)




