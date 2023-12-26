# Install required packages
# install.packages(c("glmnet", "pROC"))

# Load required libraries
library(glmnet)
library(pROC)

# Set file paths and working directory
expFile = "normalize.txt"  # Expression data file
geneFile = "score.csv"     # Gene score file
setwd("E:\\桌面\\骨关节炎单基因MR\\12.ROC")  # Set working directory

# Read expression data file
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Extract phenotype information (assuming the column names contain "_Control" or "_Case")
y = ifelse(grepl("_Control", colnames(rt)), 0, 1)

# Read gene score file
geneRT = read.csv(geneFile, header = TRUE, sep = ",", check.names = FALSE)
colnames(geneRT) = c(colnames(geneRT)[2:ncol(geneRT)], "N")
geneRT = geneRT[order(geneRT$Degree, geneRT$Betweenness, decreasing = c(TRUE, TRUE)), ]
geneRT = geneRT[1:5, ]

# Save top 10 hub genes to file
hubGene = data.frame(id = row.names(geneRT), geneRT)
hubGene = hubGene[, c("id", "Degree")]
write.table(hubGene, file = "hubGenes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Set colors for the hub genes
bioCol = rainbow(nrow(geneRT), s = 0.9, v = 0.9)

# Loop through genes and generate ROC curves
aucText = c()
rocSigGenes = c()
k = 0

for (x in row.names(geneRT)) {
  k = k + 1
  # Generate ROC curve
  roc1 = roc(y, as.numeric(rt[x,]))
  
  # Plot ROC curve
  if (k == 1) {
    pdf(file = "ROC.genes.pdf", width = 6, height = 5.5)
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "")
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", add = TRUE)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}

# Plot legend and save ROC plot
legend("bottomright", aucText, lwd = 2, bty = "n", col = bioCol[1:(ncol(rt) - 1)])
dev.off()