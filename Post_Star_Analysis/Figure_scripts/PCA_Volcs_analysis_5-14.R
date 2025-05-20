# 1. Load libraries
library(limma)
library(Glimma)

library(edgeR)
library(EnhancedVolcano)
library(dplyr)

# 2. Read counts and sample info
counts <- read.csv("./PCA_Lines_Counts.csv", row.names = 1)
colData <- read.csv("./PCA_BC_lines_Coldata.CSV", row.names = 1)


# 3. Reorder colData to match count matrix columns
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check

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


cpm_data <- cpm(dge)
keep <- rowSums(cpm_data > 2) >= (ncol(dge) / 2)  # must be >1 CPM in ≥ half the samples
dge <- dge[keep, , keep.lib.sizes = TRUE]


dge <- calcNormFactors(dge, method = "TMM")


dge$samples$lib.size            # raw total count
dge$samples$norm.factors        # scaling factor


# Run after you’ve defined dge
plotMDS(dge, labels = colData$condition, col = as.numeric(factor(colData$condition)))

# 6. Voom transformation
# 6. Voom transformation with sample quality weights
v <- voomWithQualityWeights(dge, design, plot = TRUE)


# Save voom mean-variance trend as TIFF
tiff("voom_mean_variance.tiff", width = 6, height = 5, units = "in", res = 300)

# Rerun voom to regenerate the plot
v <- voom(dge, design, plot = TRUE)

dev.off()
# 7. Fit linear model and empirical Bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)



# 8. Extract results
res <- topTable(fit, coef = "groupGroup1", number = Inf, sort.by = "logFC")
colnames(fit$coefficients)
res$symbol <- rownames(res)
res$Significant <- ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1, "Yes", "No")


# 9. Save to CSV
write.csv(res, file = "./PCA_groups_res.csv", row.names = TRUE)

#################################################################################

##################################################################################

#build XY data table
# Set a log2 fold change threshold of 1
de_status <- decideTests(fit, lfc = 1, p.value = 0.05)
summary(de_status)
# Create the volcano plot with Glimma

glimmaVolcano(fit,
              coef = "groupGroup1",
              status = de_status,
              names = rownames(fit)
)

glimmaMA(fit, coef = "groupGroup2", html = "ma_plot.html")

# How many genes pass just the FDR threshold?
sum(res$adj.P.Val < 0.05)

# How many pass logFC > 1?
sum(abs(res$logFC) > 1)

# How many pass both FDR < 0.05 & abs(logFC) > 1?
sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)

# How many pass your strict filter of FDR < 0.01 & abs(logFC) > 2?
sum(res$adj.P.Val < 0.01 & abs(res$logFC) > 2)
