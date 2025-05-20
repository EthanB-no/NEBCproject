# Load necessary libraries
library(org.Mm.eg.db)  # use org.Mm.eg.db for mouse
library(dplyr)

# Load your count matrix (rows = gene symbols, columns = samples)
counts <- read.csv("Desktop/BioInformatics/WDLSTAR_counts_matrix(in).csv", row.names = 1)

# Get all  protein-coding gene symbols
gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = rownames(counts),
                                   columns = c("SYMBOL", "GENETYPE"),
                                   keytype = "SYMBOL")

# Filter for protein-coding
protein_coding_symbols <- gene_info %>%
  filter(GENETYPE == "protein-coding") %>%
  pull(SYMBOL) %>%
  unique()

# Subset your count matrix
filtered_counts <- counts[rownames(counts) %in% protein_coding_symbols, ]

# Optional: save result
write.csv(filtered_counts, "pc_counts_5-20.csv")
