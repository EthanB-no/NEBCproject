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
    
    
    
    # 8. Extract results
    res <- topTable(fit, coef = "groupFOXA2oe", number = Inf, sort.by = "logFC")
    colnames(fit$coefficients)
    res$symbol <- rownames(res)
    res$Significant <- ifelse(res$adj.P.Val < 0.01 & abs(res$logFC) > 2, "Yes", "No")
    
    
    # 9. Save to CSV
    write.csv(res, file = "./WDL_FoxA1ShortSH_vs_NTC_limma_voom.csv", row.names = TRUE)
    
    #################################################################################
    
    ##################################################################################
    
    #build XY data table
    # Set a log2 fold change threshold of 1
    de_status <- decideTests(fit, lfc = 1, p.value = .1)
    summary(de_status)
    # Create the volcano plot with Glimma
    
    glimmaVolcano(fit,
                  coef = "groupGroup1",
                  status = de_status,
                  names = rownames(fit)
                 )
  
  glimmaMA(fit, coef = "groupFOXA2oe", html = "ma_plot.html")
  
  
 
  # Load required libraries
  library(ggplot2)
  library(ggrepel)
  
  # Assuming res_shrink is your limma result table
  
  # Define genes you want to label
  genes_of_interest <- c("Foxa2", "Tp53", "Myc")  # change these as needed
  res$gene <- rownames(res)
  # Add differential expression status
    res$DE_Status <- "Not DE"
    res$DE_Status[res$adj.P.Val < 0.01 & res$logFC > 2] <- "Over"
    res$DE_Status[res$adj.P.Val < 0.01 & res$logFC < -2] <- "Under"
    res$DE_Status <- factor(res$DE_Status, levels = c("Under", "Not DE", "Over"))
    
    # Glimma color scheme
    glimma_colors <- c("Under" = "blue", "Not DE" = "grey", "Over" = "red")
    
    # Save to JPEG
    jpeg("volcano_plot_glimma_style.jpg", width = 5, height = 5, units = "in", res = 500)
    
    ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = DE_Status), alpha = 0.6) +
      scale_color_manual(values = glimma_colors) +
      geom_text_repel(data = subset(res, gene %in% genes_of_interest),
                      aes(label = gene),
                      size = 3, max.overlaps = Inf) +
      theme_minimal(base_size = 14) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
      theme(legend.position = "top")
    
    dev.off()
    
 
