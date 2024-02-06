# setwd("/home/giang/R")
library("dplyr")
library("DESeq2")
library("purrr")
# df_list <- list(VN_0005, VN_0007, VN_0012,VN_0085,VN_0188, SRR14431113,SRR14431114,SRR14431119,SRR14431121,SRR14431122)
# salmon_count <- reduce(df_list, full_join, by = "Name")

# import data
salmon_count<- read.csv("BC_normal_count.csv",header = T)
coldata<-read.csv("BC_normal_metadata.csv")
View(coldata)
# prepare count matrix, keep salmon_count
salmon_count2=as.data.frame(salmon_count[,-1])
# row.names(salmon_count2)=make.names(salmon_count$gene_name, unique = T)
# salmon_count2<-salmon_count2 %>% dplyr::select(-gene_name)

#### TH added 10/10/2023 --- aggregate by gene names
salmon_count2 <- aggregate(salmon_count2[,2:14], by = list(salmon_count2$gene_name), mean)
rownames(salmon_count2) <- salmon_count2$Group.1
salmon_count2 <- salmon_count2[,-1]
# salmon_count2 <- log2(salmon_count2 + 1)
# boxplot(salmon_count2, las = 2, ylab = "gene expression", main = "Gene expression of BRCA case and control samples")

View(salmon_count2)
# filter low expression gene
# salmon_count2 = salmon_count2[which(rowSums(salmon_count2) > 50),]
# View(salmon_count2)

# put dataset into DEseq2
dds <- DESeqDataSetFromMatrix(countData = round(salmon_count2),
                              colData = coldata,
                              design = ~ condition)
dds
#run Deseq2
# estimating size factors: library and count nornalize
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

##TH added 10/10/2023
library(apeglm)
library("clusterProfiler")
# library(org.Hs.eg.db)
library("AnnotationDbi")
library("DESeq2")

dds <- DESeq(dds)
# extracts result table with log2 fold changes, p values and adjusted p values
res <- results(dds)
res
res <- results(dds, name= "condition_UTV_vs_normal")
res <- results(dds, contrast=c("condition","UTV","normal"))

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_UTV_vs_normal", type="apeglm")

resOrdered <- res[order(res$pvalue),]
summary(res)
res05<- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj <0.05, na.rm = T)
# plot counts

plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
# export results to csv file
write.csv(as.data.frame(resOrdered), 
          file="condition_UTV_results.csv")

####Enrichment analysis --- TH added 13/10/2023
source("../Rcode/statFuncs.R")
df <- as.data.frame(resOrdered)
df2 <- na.omit(df)
upGenes <- rownames(df2[df2$log2FoldChange > 2 & df2$pvalue < 0.05,])
downGenes <- rownames(df2[df2$log2FoldChange < -2,])

##load pathways
load("../Rdata/keggL.rda")
load("../Rdata/immuneL.rda")
load("../Rdata/hallmarkL.rda")
load("../Rdata/goL.rda")

write.csv(enrichedPW(hallmarkL, nrow(df2), upGenes, 0.1, 1)[1:5,], file = "Result/up_hallmark.csv")
write.csv(enrichedPW(immuneL, nrow(df2), upGenes, 0.1, 1)[1:5,], file = "Result/up_immune.csv")
write.csv(enrichedPW(keggL, nrow(df2), upGenes, 0.1, 1)[1:5,], file = "Result/up_kegg.csv")
write.csv(enrichedPW(goL, nrow(df2), upGenes, 0.1, 1)[1:5,], file = "Result/up_go.csv")

write.csv(enrichedPW(hallmarkL, nrow(df2), downGenes, 0.1, 1)[1:5,], file = "Result/down_hallmark.csv")
write.csv(enrichedPW(immuneL, nrow(df2), downGenes, 0.1, 1)[1:5,], file = "Result/down_immune.csv")
write.csv(enrichedPW(keggL, nrow(df2), downGenes, 0.1, 1)[1:5,], file = "Result/down_kegg.csv")
write.csv(enrichedPW(goL, nrow(df2), downGenes, 0.1, 1)[1:5,], file = "Result/down_go.csv")

################################
################################
#transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dds)

meanSdPlot(assay(ntd))

# clustering and visualization
vsdata<- vst(dds, blind = F, fitType = 'local')
dim(vsdata)
plotPCA(vsdata, intgroup = "condition")
# plot disperation
plotDispEsts(dds)


# heatmap
sig.df<- as.data.frame(res05)
View(sig.df)


# filtering
sig.df<- sig.df[(sig.df$baseMean>150)&(abs(sig.df$log2FoldChange)>1.5),]

mat<- counts(dds, normalized =T)[rownames(sig.df),]
mat.z<- t(apply(mat,1,scale))
mat.z
colnames(mat.z)<- coldata$sample
# Define a custom color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap(mat.z, cluster_rows= T, cluster_columns=T, column_labels=colnames(mat.z),
        col = my_palette)

# volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# GO
genes_set<-rownames(sig.df)

GO_results <- enrichGO(gene = gene_set, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 15))


