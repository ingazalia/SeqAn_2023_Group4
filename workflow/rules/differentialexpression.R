#RNA SEQ
#chooseCRANmirror()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")#,lib="results/diffexpression")
if(!require(DESeq2)){
  BiocManager::install('DESeq2')#, repos = "http://cran.us.r-project.org")#,lib="results/diffexpression")
}
if(!require(pheatmap)){
  #install.packages("pheatmap")
  install.packages("pheatmap", repos = "http://cran.us.r-project.org")#,lib="results/diffexpression")
}

library(DESeq2)
library(pheatmap)

source('workflow/scripts/differentialexpression_function.R')


# import data
rnaseqdata <- read.table(snakemake@input[[1]], sep = "\t",header = TRUE,row.names = 1,check.names=F)
sample_sheet <- read.table(snakemake@input[[2]],header = TRUE,row.names = 1)

print('Dim of FeatureCount data')
dim(rnaseqdata)
# 5527 different genes and 5 samples
print('Any NAs?')
any(is.na(rnaseqdata))

# Groups
group_A <- rnaseqdata[,rownames(sample_sheet)[which(sample_sheet$group=='A')]]
colnames(group_A) <- paste0(colnames(group_A),'_A')
group_B <- rnaseqdata[,rownames(sample_sheet)[which(sample_sheet$group=='B')]]
colnames(group_B) <- paste0(colnames(group_B),'_B')
                            
# new order of patient
rnaseqdata <- data.frame(group_A,group_B)

# First task: Comparing gene expression between AD and old samples

# DIFFERENTIAL EXPRESSION
res <- diffexpression(group_A,group_B,alpha=0.05,tidy=FALSE)
res <- res[!is.na(res$log2FoldChange),]
summary(res)
dds <- diffexpression(group_A,group_B,alpha=0.05,tidy=FALSE,result=FALSE)

resOrdered <- res[order(res$padj),]
#Yekutieli control 

#multiple test control
pvalues <- p.adjust(resOrdered$pvalue, method = "BY")
print('How many genes are significant diff. expressed?')
sum(pvalues<0.05,na.rm = T)


plotdf<- data.frame(genes=rownames(resOrdered),log2FoldChange=resOrdered$log2FoldChange, padj= resOrdered$padj, Yekutielipval=pvalues,pval=resOrdered$pvalue)
#plotdf <- plotdf[!is.na(plotdf$log2FoldChange),]
plotdf$gene_type <- 'not significant'
plotdf$gene_type[plotdf$log2FoldChange > 0 & plotdf$padj < 0.05] <- "upregulated"
plotdf$gene_type[plotdf$log2FoldChange < 0 & plotdf$padj < 0.05] <- "downregulated"
print("Table # down/upregulated genes")
print(table(plotdf$gene_type))

write.table(plotdf,snakemake@output[[1]],quote = FALSE, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)

#up and downregulated genes, ordered by p_adjust
upregulated <- plotdf[plotdf$gene_type=="upregulated",]
upregulated <- upregulated[order(upregulated$padj),]
downregulated <- plotdf[plotdf$gene_type=="downregulated",]
downregulated <- downregulated[order(downregulated$padj),]

# Save up-/downregulated count matrix and up-/downregulated genes 
write.table(as.factor(unlist(upregulated$genes)), snakemake@output[[2]],quote = FALSE, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)

write.table(upregulated,snakemake@output[[3]] , quote=FALSE, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)


write.table(downregulated$genes, snakemake@output[[4]],quote = FALSE, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(downregulated, snakemake@output[[5]], quote = FALSE, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)


################################################################################
## VISUALIZAION
pdf("results/diffexpression/plots_combined.pdf")
#setwd("/")
# Plot counts for a single gene -> smallest pvalue
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# MA plot
png(file ="results/diffexpression/MAplot.png" , height = 400,width=400)
plotMA(res, ylim=c(-3,3),colSig='red',colLine='red',alpha=0.05,cex=0.4,xlab='Mean expression',ylab='log2 (Group B/Group A)')
dev.off()
# description: MA plot of differential results (AD vs old)

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
plotMA(resNorm, ylim=c(-1.5,1.5),main='Normal')



#Heatmap of the count matrix + cluster
upregulated <- upregulated[order(upregulated$log2FoldChange,decreasing = TRUE),]
downregulated <- downregulated[order(downregulated$log2FoldChange),]
ntd <- assay(normTransform(dds))
ntd <- assay(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:60]
png(file = "results/diffexpression/pheatmap_up.png", height = 600,width=400)
pheatmap(ntd[which(rownames(ntd) %in% upregulated$genes[1:10]),],cluster_cols = FALSE)
dev.off()

png(file = "results/diffexpression/pheatmap_down.png", height = 600,width=400)
pheatmap(ntd[which(rownames(ntd) %in% downregulated$genes[1:15]),],cluster_cols = FALSE)
dev.off()

##Principal component plot of the samples
#Related to the distance matrix is the PCA plot, which shows the samples in the
#2D plane spanned by their first two principal components. This type of plot is
#useful for visualizing the overall effect of experimental covariates and batch effects.
vsd <- vst(dds, blind=FALSE)
vsd$condition <- factor(c(rep("A",ncol(group_A)),rep("B",ncol(group_B))))
levels(vsd$condition) <- c('A','B')
png(file = "results/diffexpression/PCARNAseq.png", height = 400,width=600)
plotPCA(vsd, intgroup=c("condition"))
dev.off()
dev.off()
# description: PCA plot A vs B 

