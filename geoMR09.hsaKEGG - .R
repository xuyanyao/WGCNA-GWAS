# Install required packages
# install.packages(c("colorspace", "stringi", "ggplot2", "circlize", "RColorBrewer"))

# Install Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
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
library(ComplexHeatmap)

# Set parameters
pvalueFilter = 0.05     # p-value filter
qvalueFilter = 1        # q-value filter
inputFile = "interGenes.txt"     # Input gene list file

# Set color selection based on q-value
colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

# Set working directory
setwd("E:\\桌面\\177geoMR\\09.KEGG")

# Read gene list
rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)

# Convert gene symbols to entrez IDs
genes = unique(as.vector(rt[, 1]))
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs = as.character(entrezIDs)
rt = data.frame(genes, entrezID = entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]

# Perform KEGG pathway enrichment analysis
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
kk@result$Description = gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG = as.data.frame(kk)
KEGG$geneID = as.character(sapply(KEGG$geneID, function(x) paste(rt$genes[match(strsplit(x, "/")[[1]], as.character(rt$entrezID))], collapse = "/")))
KEGG = KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter),]

# Write KEGG results to file
write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Set the number of pathways to display
showNum = 30     # Display the top 30 pathways
if (nrow(KEGG) < showNum) {
  showNum = nrow(KEGG)
}

# Generate bar plot
pdf(file = "barplot.pdf", width = 8.5, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 130, color = colorSel)
dev.off()

# Generate bubble plot
pdf(file = "bubble.pdf", width = 8.5, height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format = 130, color = colorSel)
dev.off()
