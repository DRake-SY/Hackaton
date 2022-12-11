## Installation des d?pendances ----------
require(devtools)
list.of.packages <- c("DESeq2", "apeglm", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(new.packages)
}

library(DESeq2)
library(apeglm)
library(dplyr)

## R?cup?ration des donn?es ------------
fic <- read.table("result.counts")
#formatage des donn?es
cnts <- as.data.frame(apply(fic[2:nrow(fic),7:length(fic)], 2, as.numeric))
colnames(cnts) <- fic[1,7:length(fic)]
row.names(cnts) <- fic[2:nrow(fic),1]
dim(cnts)

## Analyse diff?rentielle primaire -------

# cond <- factor(c(1,1,1,2,2,2,2,2)) # diff?rents goupes : 1 = mut?, 2 = non-mut?
# dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
# 
# # standard analysis
# dds <- DESeq(dds)
# res <- results(dds)
# summary(dds)
# summary(res)
# head(res)

# # moderated log2 fold changes
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# 
# # an alternate analysis: likelihood ratio test
# ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
# resLRT <- results(ddsLRT)
# 
# par(mfrow=c(3,3))
# hist(res$pvalue)
# hist(resLFC$pvalue)
# hist(resLRT$pvalue)
# hist(res$padj)
# hist(resLFC$padj)
# hist(resLRT$padj)
# hist(res$log2FoldChange)
# hist(resLFC$log2FoldChange)
# hist(resLRT$log2FoldChange)
# 
# dds[res$pvalue[1-is.na(res$pvalue)]]
# 
# dds[resLFC$pvalue<0.01]
# 
# dds[resLRT$pvalue<0.01]

#les r?sultats sont inexploitables
#on va filtrer les donn?es afin de pouvoir refaire l'analyse

## Filtrage des g?nes -------

# test de diff?rents filtres

#on enl?ve tous les g?nes qui au moins une fois ne sont pas exprim?s
cnts_filtered_nozero <- subset(cnts, cnts[,1] > 0 & cnts[,2] > 0 & cnts[,3] > 0 & cnts[,4] > 0 & cnts[,5] > 0 & cnts[,6] > 0 & cnts[,7] > 0 & cnts[,8] > 0)
#seuillage, enlever les g?nes qui en moyenne ne sont pas exprim?s plus de 5 fois
cnts_filtered_mean5 = cnts[which(rowMeans(cnts)>5),]
#on enl?ve les g?nes non exprim?s
#cnts_filtered_mean0 = cnts[-which(rowMeans(cnts)==0),]

dim(cnts_filtered_nozero)
# 15839 g?nes gard?s
#summary(cnts_filtered_nozero)
dim(cnts_filtered_mean5)
# 19210 g?nes gard?s
#summary(cnts_filtered_mean5)
#dim(cnts_filtered_mean0)
# 30029 g?nes gard?s
#summary(cnts_filtered_mean0)

## Analyse diff?rentielle sur table filtr?e --------

# Construction d'un objet DEseq
cond <- factor(c(1,1,1,2,2,2,2,2)) # diff?rents goupes : 1 = mut?, 2 = non-mut?
dds_nozero <- DESeqDataSetFromMatrix(cnts_filtered_nozero, DataFrame(cond), ~ cond)
dds_mean5 <- DESeqDataSetFromMatrix(cnts_filtered_mean5, DataFrame(cond), ~ cond)
#dds_mean0 <- DESeqDataSetFromMatrix(cnts_filtered_mean0, DataFrame(cond), ~ cond)

# standard analysis
dds_nozero = DESeq(dds_nozero)
dds_mean5 = DESeq(dds_mean5)
#dds_mean0 = DESeq(dds_mean0)
res_nozero = results(dds_nozero)
res_mean5 = results(dds_mean5)
#res_mean0 = results(dds_mean0)

par(mfrow=c(1,3))
hist(res_nozero$pvalue, main="nozero")
#hist(res_mean0$pvalue, main="mean0")
hist(res_mean5$pvalue, main="mean5")
#le seuillage ? 0 est insuffisant, les p-values sont trop h?t?rog?nes
#on arr?te ici l'analyse pour ce jeu de donn?es

# moderated log2 fold changes
resLFC_mean5 <- lfcShrink(dds_mean5, coef=2, type="apeglm")
resLFC_nozero <- lfcShrink(dds_nozero, coef=2, type="apeglm")

# an alternate analysis: likelihood ratio test
ddsLRT_mean5 <- DESeq(dds_mean5, test="LRT", reduced= ~ 1)
resLRT_mean5 <- results(ddsLRT_mean5)
ddsLRT_nozero <- DESeq(dds_nozero, test="LRT", reduced= ~ 1)
resLRT_nozero <- results(ddsLRT_nozero)

pdf("p_nozero.pdf")
par(mfrow=c(1,3))
hist(res_nozero$pvalue, main='res pvalues')
hist(resLFC_nozero$pvalue, main='LFC pvalues')
hist(resLRT_nozero$pvalue, main='LRT pvalues')
dev.off()
pdf("p_adj_nozero.pdf")
par(mfrow=c(1,3))
hist(res_nozero$padj, main='res padj')
hist(resLFC_nozero$padj, main='LFC padj')
hist(resLRT_nozero$padj, main='LRT padj')
dev.off()
pdf("LogFC_no_zero.pdf")
par(mfrow=c(1,3))
hist(res_nozero$log2FoldChange, main='res logFC')
hist(resLFC_nozero$log2FoldChange, main='LFC logFC')
hist(resLRT_nozero$log2FoldChange, main='LRT logFC')
dev.off()
pdf("p_mean5.pdf")
par(mfrow=c(1,3))
hist(res_mean5$pvalue, main='res pvalues')
hist(resLFC_mean5$pvalue, main='LFC pvalues')
hist(resLRT_mean5$pvalue, main='LRT pvalues')
dev.off()
pdf("p_adj_mean5.pdf")
par(mfrow=c(1,3))
hist(res_mean5$padj, main='res padj')
hist(resLFC_mean5$padj, main='LFC padj')
hist(resLRT_mean5$padj, main='LRT padj')
dev.off()
pdf("LogFC_mean5.pdf")
par(mfrow=c(1,3))
hist(res_mean5$log2FoldChange, main='res logFC')
hist(resLFC_mean5$log2FoldChange, main='LFC logFC')
hist(resLRT_mean5$log2FoldChange, main='LRT logFC')
dev.off()

# dds[res$pvalue[1-is.na(res$pvalue)]]
# 
# dds[resLFC$pvalue<0.01]
# 
# dds[resLRT$pvalue<0.01]