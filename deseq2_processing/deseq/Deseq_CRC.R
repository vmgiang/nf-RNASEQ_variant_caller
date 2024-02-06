setwd("/home/giang/R/DEG/DEG_CRC/")
library("dplyr")
library("DESeq2")

# import read count
salmon_count<- read.csv("CRC_salmon_count.csv",header = T)
View(salmon_count)
# inport meta data
coldata<-read.csv("CRC_meta_data.csv")
View(coldata)
# prepare count matrix, keep salmon_count
salmon_count2=as.data.frame(salmon_count[,-1])
# row.names(salmon_count2)=make.names(salmon_count$gene_name, unique = T)
# salmon_count2<-salmon_count2 %>% dplyr::select(-gene_name)

#### TH added 10/10/2023 --- aggregate by gene names
salmon_count2 <- aggregate(salmon_count2[,2:20], by = list(salmon_count2$gene_name), mean)
rownames(salmon_count2) <- salmon_count2$Group.1
salmon_count2 <- salmon_count2[,-1]
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
res <- results(dds, name= "condition_normal_vs_CRC")
res <- results(dds, contrast=c("condition","normal","CRC"))

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_normal_vs_CRC", type="apeglm")

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
          file="condition_CRC_results.csv")

####Enrichment analysis --- TH added 13/10/2023
source("../Rcode/statFuncs.R")

df <- as.data.frame(resOrdered)
df2 <- na.omit(df)
upGenes <- rownames(df2[df2$log2FoldChange > 12 & df2$pvalue < 0.05,])
downGenes <- rownames(df2[df2$log2FoldChange < -12,])
write.csv(upGenes, file = "upGenes.csv",sep = ",")
write.csv(downGenes,file = "downGenes.csv",sep = ",")
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
sig.df<- na.omit(as.data.frame(res05))


# filtering
sig.df<- sig.df[(sig.df$baseMean>200)&(abs(sig.df$log2FoldChange)>8),]

mat<- counts(dds, normalized =T)[rownames(sig.df),]
mat.z<- t(apply(mat,1,scale))
mat.z
colnames(mat.z)<- coldata$sample
# Define a custom color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap(mat.z, cluster_rows= T, cluster_columns=T, column_labels=colnames(mat.z),
        col = my_palette)


library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Add threshold lines
ggplot(data = df2, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point() 
# Biostatsquid theme
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
theme_minimal()

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"
# if log2FoldChange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$log2FoldChange > 12 & df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -12 & df$pvalue < 0.05] <- "DOWN"

head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])
head(df[order(df$padj) & df$diffexpressed == 'UP', ])
# Create a new column "delabel" to de, that will contain the name of the top 5 differentially expressed genes (NA in case they are not)
df$gene_symbol<- rownames(df)
df$delabel <- ifelse((df$gene_symbol %in% head(df[order(df$padj) & df$diffexpressed == 'UP', "gene_symbol"], 5)) |
                       (df$gene_symbol %in% head(df[order(df$padj) & df$diffexpressed == 'DOWN', "gene_symbol"], 5)),
                     df$gene_symbol, NA)

ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 70), xlim = c(-30, 30)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-30, 30, 5)) + # to customise the breaks in the x axis
  ggtitle('Dysregualted genes in Colorectal cancer samples vs Normal samples') + # Plot title 
  geom_text_repel(max.overlaps = Inf) 

