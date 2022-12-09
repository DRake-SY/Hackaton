library(DESeq2)
library(apeglm)
library(dplyr)
# see vignette for suggestions on generating
# count tables from RNA-Seq data
# modifier l'input lorsque l'on aura les tables de comptage
fic <- read.table("result.counts")
cnts <- as.data.frame(apply(fic[2:nrow(fic),7:length(fic)], 2, as.numeric))
colnames(cnts) <- fic[1,7:length(fic)]
row.names(cnts) <- fic[2:nrow(fic),1]
dim(cnts)
fil <- 10 ; 
# filtrage
# cnts <- subset(cnts, cnts[,1] > fil | cnts[,2] > fil | cnts[,3] > fil | cnts[,4] > fil | cnts[,5] > fil | cnts[,6] > fil | cnts[,7] > fil | cnts[,8] > fil); : on enleve la ligne si elle est nulle
# cnts <- subset(cnts, cnts[,1] > fil & cnts[,2] > fil & cnts[,3] > fil & cnts[,4] > fil & cnts[,5] > fil & cnts[,6] > fil & cnts[,7] > fil & cnts[,8] > fil); : on enlève la ligne si elle contient un zero
cnts <- subset(cnts, mean(cnts[,1:3]) > fil & mean(cnts[,4:8]) > fil);
# cnts <- subset(cnts, mean(cnts[,1:3]) > fil | mean(cnts[,4:8]) > fil);
dim(cnts)
summary(cnts)
cond <- factor(c(1,1,1,2,2,2,2,2)) # différents goupes : 1 = muté, 2 = non-muté
summary(cond)

# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
head(dds)

# standard analysis
dds <- DESeq(dds)
res <- results(dds)
summary(dds)
summary(res)
head(res)

# moderated log2 fold changes
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

# an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)

par(mfrow=c(3,3))
hist(res$pvalue)
hist(resLFC$pvalue)
hist(resLRT$pvalue)
hist(res$padj)
hist(resLFC$padj)
hist(resLRT$padj)
hist(res$log2FoldChange)
hist(resLFC$log2FoldChange)
hist(resLRT$log2FoldChange)

dds[res$pvalue[1-is.na(res$pvalue)]]

dds[resLFC$pvalue<0.01]

dds[resLRT$pvalue<0.01]

# PCA sur les données filtrées