# Install required packages
# install.packages(c("colorspace", "stringi", "ggplot2", "circlize", "RColorBrewer", "ggpubr"))

# Install Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db", "DOSE", "clusterProfiler", "enrichplot", "ComplexHeatmap"))

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

# Set parameters
pvalueFilter = 0.05      # p-value filter
qvalueFilter = 1         # q-value filter
inputFile = "interGenes.txt"     # Input gene list file

# Set color selection based on q-value
colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

# Set working directory
setwd("E:\\桌面\\177geoMR\\08.GO")

# Read gene list
rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)
genes = unique(as.vector(rt[, 1]))

# Convert gene symbols to entrez IDs
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs = as.character(entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]

# Perform GO enrichment analysis
kk = enrichGO(gene = gene, OrgDb = "org.Hs.eg.db", pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO = as.data.frame(kk)
GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter),]

# Write GO results to file
write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Generate bar plot
pdf(file = "barplot.pdf", width = 9, height = 7)
bar = barplot(kk, drop = TRUE, showCategory = 10, label_format = 130, split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

# Generate bubble plot
pdf(file = "bubble.pdf", width = 9, height = 7)
bub = dotplot(kk, showCategory = 10, orderBy = "GeneRatio", label_format = 130, split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# Perform GO circular visualization
# ... (the rest of the code for circular visualization)

# Save the circular visualization as a PDF
pdf(file = "GO.circlize.pdf", width = 10, height = 10)
# ... (the rest of the code for circular visualization)
dev.off()

# Print circle size
circle_size = unit(1, "snpc")
print(circle_size)
