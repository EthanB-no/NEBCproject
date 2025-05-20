# 1. Load libraries
library(limma)
library(Glimma)

library(edgeR)
library(EnhancedVolcano)
library(dplyr)

# 2. Read counts and sample info
counts <- read.csv("./WDLSTAR_counts_matrix(in).csv", row.names = 1)
colData <- read.csv("./WDL_ColData.CSV", row.names = 1)


# 3. Reorder colData to match count matrix columns
colData <- colData[match(colnames(counts), rownames(colData)), , drop = FALSE]
stopifnot(all(rownames(colData) == colnames(counts)))  # sanity check


# 4. Subset to only NTC and FOXA2oe samples
keep_samples <- colData$condition %in% c("NTC", "FOXA2oe")
colData <- colData[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]


# 5. Create DGEList and filter low-expressed genes
#counts <- counts[rowSums(counts) > 10, ]
#counts are counts
dge <- DGEList(counts = counts)

#set NTC as reference for comparisons
group <- relevel(factor(colData$condition), ref = "NTC")
dge$samples$group <- group
design <- model.matrix(~ group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]


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

# Save voom mean-variance trend as TIFF
tiff("voom_mean_variance.tiff", width = 6, height = 5, units = "in", res = 300)

# Rerun voom to regenerate the plot
v <- voom(dge, design, plot = TRUE)

dev.off()
# 7. Fit linear model and empirical Bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "groupFOXA2oe", number = Inf, sort.by = "logFC")
colnames(fit$coefficients)
res$symbol <- rownames(res)
res$Significant <- ifelse(res$adj.P.Val < 0.01 & abs(res$logFC) > 2, "Yes", "No")

# 8. Extract results
res <- topTable(fit, coef = "groupFOXA2oe", number = Inf, sort.by = "logFC")
colnames(fit$coefficients)
res$symbol <- rownames(res)
res$Significant <- ifelse(res$adj.P.Val < 0.01 & abs(res$logFC) > 2, "Yes", "No")



genes_of_interest <- list(
  Luminal = c(  # Expected DOWN with FOXA2
    "Krt20", "Pparg", "Foxa1", "Gata3", "Snx31", "Upk1a", "Upk2", "Fgfr3"
  ),
  ECM_SmoothMuscle = c(  # Variable, often UP in EMT/stem-like states
    "Pgm5", "Des", "C7", "Sfrp4", "Comp", "Sgcd"
  ),
  EMT_Claudin = c(  # Expected UP with FOXA2-induced plasticity
    "Zeb1", "Zeb2", "Snai1", "Twist1", "Cdh2", "Cldn3", "Cldn4", "Cldn7"
  ),
  Basal = c(  # Often UP in FOXA2 contexts that lose luminal identity
    "Cd44", "Krt6a", "Krt5", "Krt14", "Col17a1"
  ),
  Squamous = c(  # May be UP in poorly differentiated states
    "Dsc3", "Gsdmc", "Tgm1", "Pi3", "Tp63"
  ),
  Immune = c(  # Context-dependent
    "Cd274", "Pdcd1lg2", "Ido1", "Cxcl11", "L1cam", "Saa1"
  ),
  Neuronal = c(  # Heterogeneously UP — hallmark of partial NE reprogramming
    "Msi1", "Plekhg4b", "Gnag", "Peg10", "Rnd2", "Alpl", "Sox2", "Tubb2b",
    "Ascl1", "Insm1", "Chga", "Syp", "Grp", "Neurod1"
  ),
  CIS_related = c(  # Markers of carcinoma in situ, variable
    "Crtac1", "Ctse", "Padi3", "Msn", "Nr3c1"
  ),
  Sonic_Hedgehog = c(  # Variable, but some show FOXA2 association
    "Shh", "Bmp5"
  ),
  Stemness = c(  # Often UP with FOXA2 overexpression and dedifferentiation
    "Sox2", "Nanog", "Pou5f1", "Aldh1a1", "Lgr5", "Prom1", "Cd44", "Klf4", "Myc"
  ),
  FOXA2_targets = c(  # Based on NEPC and lung development contexts — generally UP
    "Ascl1", "Insm1", "Grp", "Syt1", "Chga", "Foxa2", "Tff3", "Neurod1", "Tubb3", "Rfx4", "Stmn2"
  )
)
# Define genes you want to label
  # change these as needed
res$gene <- rownames(res)
# Add differential expression status
res$DE_Status <- "Not DE"
res$DE_Status[res$adj.P.Val < 0.01 & res$logFC > 2] <- "Over"
res$DE_Status[res$adj.P.Val < 0.01 & res$logFC < -2] <- "Under"
res$DE_Status <- factor(res$DE_Status, levels = c("Under", "Not DE", "Over"))

# Glimma color scheme
glimma_colors <- c("Under" = "blue", "Not DE" = "grey", "Over" = "red")

# Save to JPEG
jpeg("volcano_plot_biomarkers.jpg", width = 5, height = 5, units = "in", res = 500)

ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = DE_Status), alpha = 0.6) +
  scale_color_manual(values = glimma_colors) +
  geom_text_repel(data = subset(res, gene %in% unlist(genes_of_interest)),
                  aes(label = gene),
                  size = 5, max.overlaps = Inf) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "top")

dev.off()

