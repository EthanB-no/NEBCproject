library(limma)
library(Glimma)
library(edgeR)
library(EnhancedVolcano)
library(dplyr)

#counts and colData
counts <- read.csv("./Data/pc_counts_5-20.csv", row.names = 1)
colData <- read.csv("./Data/WDL_ColData.CSV", row.names = 1)


# 3. Reorder colData to match count matrix columns
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check

# 4. Subset to only NTC and FOXA2oe samples
keep_samples <- colData$condition %in% c("NTC", "FOXA2oe")
colData <- colData[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]


dge <- DGEList(counts = counts)

#set NTC as reference for comparisons -> filter by expression plus dge steps
group <- relevel(factor(colData$condition), ref = "NTC")
dge$samples$group <- group
design <- model.matrix(~ group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = TRUE]


#normalize by trimmed mean of m values
dge <- calcNormFactors(dge, method = "TMM")

#some sanity checks here
dge$samples$lib.size            # raw total count
dge$samples$norm.factors        # scaling factor


# Run after you’ve defined dge (duh)
plotMDS(dge, labels = colData$condition, col = as.numeric(factor(colData$condition)))

# 6. Voom transformation
v <- voom(dge, design, plot = TRUE)


# Save voom mean-variance trend as TIFF NO NEED TO RUN UNLESS YOU WANT TO
tiff("voom_mean_variance.tiff", width = 6, height = 5, units = "in", res = 300)

# Rerun voom to regenerate the plot
v <- voom(dge, design, plot = TRUE)
dev.off()


# 7. Fit linear model and empirical Bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)


# 8. Extract results
res <- topTable(fit, coef = "groupFOXA2oe", number = Inf, sort.by = "logFC")
colnames(fit$coefficients)
res$symbol <- rownames(res) 
res$Significant <- ifelse(res$adj.P.Val < 0.01 & abs(res$logFC) > 2, "Yes", "No")


# 9. Save to CSV + annotate
library(org.Mm.eg.db)
library(AnnotationDbi)

# Map gene symbols to full gene names (e.g., full descriptions)
gene_info <- select(org.Mm.eg.db,
                    keys = res$symbol,
                    columns = c("SYMBOL", "GENENAME"),
                    keytype = "SYMBOL")

# Join the result with your res table
res_annotated <- left_join(res, gene_info, by = c("symbol" = "SYMBOL"))
res_annotated <- res_annotated %>%
  dplyr::select(symbol, everything())


# 9. Save to CSV
write.csv(res_annotated, file = "./5-21_Foxa2_res.csv", row.names = FALSE)

#################################################################################
#plotting!
##################################################################################

#build XY data table
# Set a log2 fold change threshold of 1
de_status <- decideTests(fit, lfc = 1.5, p.value = 0.01)
summary(de_status)


# Create the volcano plot with Glimma

glimmaVolcano(fit,
              coef = "groupFOXA2oe",
              status = de_status,
              names = rownames(fit)
)

glimmaMA(fit, coef = "groupGroup2", html = "ma_plot.html")

# How many genes pass just the p.adj threshold?
sum(res$adj.P.Val < 0.05)

# How many pass logFC > 1?
sum(abs(res$logFC) > 1)

# How many pass both FDR < 0.05 & abs(logFC) > 1?
sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)

# How many pass your strict filter of FDR < 0.01 & abs(logFC) > 2?
sum(res$adj.P.Val < 0.01 & abs(res$logFC) > 2)
