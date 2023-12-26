# Load libraries
library(limma)
library(pheatmap)
library(ggplot2)

# Set input file paths and parameters
inputFile <- "geneMatrix.txt"
conFile <- "s1.txt"
treatFile <- "s2.txt"
logFCfilter <- 2
adj.P.Val.Filter <- 0.05


# Read gene expression data
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt <- avereps(data)

# Log transformation
qx <- as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)
if (LogC) {
  rt[rt < 0] <- 0
  rt <- log2(rt + 1)
}
data <- normalizeBetweenArrays(rt)

# Read sample information
sample1 <- read.table(conFile, header = FALSE, sep = "\t", check.names = FALSE)
sample2 <- read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
sampleName1 <- gsub("^ | $", "", as.vector(sample1[, 1]))
sampleName2 <- gsub("^ | $", "", as.vector(sample2[, 1]))
conData <- data[, sampleName1]
treatData <- data[, sampleName2]
data <- cbind(conData, treatData)
conNum <- ncol(conData)
treatNum <- ncol(treatData)

# Design matrix and linear modeling
Type <- c(rep("con", conNum), rep("treat", treatNum))
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("con", "treat")
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(treat - con, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Top differentially expressed genes
allDiff <- topTable(fit2, adjust = 'BY', number = 200000)
allDiffOut <- rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "all.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Save normalized data
Type <- c(rep("Control", conNum), rep("Treat", treatNum))
outData <- rbind(id = paste0(colnames(data), "_", Type), data)
write.table(outData, file = "normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Extract significantly differentially expressed genes
diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = "diff.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Generate heatmap
geneNum <- 50
diffUp <- diffSig[diffSig$logFC > 0, ]
diffDown <- diffSig[diffSig$logFC < 0, ]
geneUp <- row.names(diffUp)
geneDown <- row.names(diffDown)
geneUp <- if (nrow(diffUp) > geneNum) geneUp[1:geneNum] else geneUp
geneDown <- if (nrow(diffDown) > geneNum) geneDown[1:geneNum] else geneDown
hmExp <- data[c(geneUp, geneDown), ]

Type <- c(rep("Control", conNum), rep("Treat", treatNum))
names(Type) <- colnames(data)
Type <- as.data.frame(Type)

pdf(file = "heatmap.pdf", width = 8, height = 6.5)
pheatmap(hmExp,
         annotation_col = Type,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 7,
         fontsize_row = 5,
         fontsize_col = 7)
dev.off()

# Volcano plot
allDiff$logFC[allDiff$logFC > 20] <- 20
allDiff$logFC[allDiff$logFC < -20] <- -20
Significant <- ifelse((allDiff$adj.P.Val < adj.P.Val.Filter & abs(allDiff$logFC) > logFCfilter),
                      ifelse(allDiff$logFC > logFCfilter, "Up", "Down"), "Not")

p <- ggplot(allDiff, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(col = Significant)) +
  scale_color_manual(values = c("green", "grey", "red")) +
  labs(title = " ") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p <- p + theme_bw()

pdf(file = "vol.pdf", width = 5.5, height = 5)
print(p)
dev.off()
