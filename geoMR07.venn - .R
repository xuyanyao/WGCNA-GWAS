# Install required package
# install.packages("VennDiagram")

# Load library for Venn diagram
library(VennDiagram)

# Set file paths and working directory
diffFile <- "diff.txt"
wgcnaFile <- "module_blue.txt"
setwd("E:\\桌面\\177geoMR\\07.venn")

# Create a list to store gene sets
geneList <- list()

# Read differential expression file
rt <- read.table(diffFile, header = TRUE, sep = "\t", check.names = FALSE)
geneNames <- as.vector(rt[, 1])
geneNames <- gsub("^ | $", "", geneNames)
uniqGene <- unique(geneNames)
geneList[["DEG"]] <- uniqGene

# Read WGCNA module file
rt <- read.table(wgcnaFile, header = FALSE, sep = "\t", check.names = FALSE)
geneNames <- as.vector(rt[, 1])
geneNames <- gsub("^ | $", "", geneNames)
uniqGene <- unique(geneNames)
geneList[["WGCNA"]] <- uniqGene

# Generate Venn diagram
color <- c("#FF0000FF", "#0000CCFF")
venn.plot <- venn.diagram(geneList, filename = NULL, fill = color, scaled = FALSE, cat.pos = c(-1, 1), cat.col = color, cat.cex = 1.2)
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# Write intersection genes to a file
interGenes <- Reduce(intersect, geneList)
write.table(interGenes, file = "interGenes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
