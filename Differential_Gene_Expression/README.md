
## RNAseq Differential Expression Analysis of _P. lutea_ adult samples from mangrove and reef sites ######

This script is based on the work of [Erin Chille](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/4-Differential-Gene-Expression-Analysis/pHTreatment_RNAseqDE.Rmd), with some modifications and adjustments.

### Set up workspace in R Studio

```{r}
# Load libraries
library("genefilter")           #BiocManager::install("genefilter") 
library("DESeq2")               #BiocManager::install("DESeq2")
library("factoextra")           #install.packages("factoextra")
library("NbClust")              #install.packages("NbClust")
library("ComplexHeatmap")       #BiocManager::install("ComplexHeatmap")
library("tidyverse")            
library("RColorBrewer")
library("ggplot2")              
library("goseq")                #BiocManager::install("goseq")
library("gridExtra")            #install.packages("gridExtra")
library("VennDiagram")          #install.packages("VennDiagram")
library("patchwork")            #install.packages("patchwork")
library("dplyr")

#load treatment information
treatmentinfo <- read.csv("RNAseq_data.csv", header = TRUE, sep = ";")
str(treatmentinfo)
head(treatmentinfo)

#load gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id")) #the gene_count_matrix.csv is the output of Stringtie
head(gcount)
```
### Construct DESeq2 dataset
#### Pre-filter gene counts
Pre-filtering our dataset to reduce the memory size dataframe, increase the speed of the transformation and improve sensitivity of statistical analysis by removing low-coverage counts. 

```
#Set filter values for PoverA: smallest sample size per treat. is 4, so 4/8 (8 samples) is 0.5
#This means that 4 out of 8 (0.5) samples need to have counts over 10. So P=50 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.5,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#Merge the age and depth columns into a new column, group. Set group as a factor.
treatmentinfo$location <- factor(treatmentinfo$location, levels = c("reef","mangrove"))

#Create a DESeqDataSet design from gene count matrix and labels. 
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                  colData = treatmentinfo,
                                  design = ~location)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
#Log-transform the data using a variance stabilizing transforamtion (vst) for visualization purposes. This transformation deals with the sampling variability of low counts by #calculating within-group variability.  

#To use vst we first need to calculate the size factors of the samples, that is an estimate of how many reads each sample has compared to the others. 
SF.gdds <- estimateSizeFactors(gdds) #size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors. In this case size factors are all less than 4, so vst can be used.

gvst <- vst(gdds, blind=FALSE) 

### Principal component plot of samples

gPCAdata <- plotPCA(gvst, intgroup = c("location"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
PCA <- ggplot(gPCAdata, aes(PC1, PC2, color=location)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(labels = c("reef","mangrove"), values = c("reef"="deepskyblue", "mangrove"="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA #set title attributes`

ggsave(file = "PCA_preDE_vst.png", PCA)
```

### Differential Gene Expression Analysis 

#### Run DE analysis - Adult samples
```{r}
#DESEq2 internally applies the median of ratios method for normalization.
DEGs <- DESeq(gdds) #run differential expression test by group using the Wald model

#Explore significant p-values for meso and shallow adults
DEGs.results <- results(DEGs, contrast= c("location","mangrove","reef"))
write.csv(DEGs.results, "DEGs_mangrove_vs_reef.csv")
head(DEGs.results)
sum(DEGs.results$padj < 0.05, na.rm=TRUE)

library("ggpubr")
DEGs.MA <- ggmaplot(
  DEGs.results,
  fdr = 0.05, size = 0.4,
  detection_call = NULL,
  genenames = NULL,
  alpha = 1,
  top = 0,
  palette = c("orange", "deepskyblue", "gray"),
  label.rectangle = FALSE,
  select.top.method = c("padj"),
  xlab = "Log2 mean expression",
  ylab = "Log2 fold change",
  ggtheme = theme_classic(),
 )

pdf(file="DEGs.MA.pdf")
DEGs.MA #view plot
dev.off()

# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DEGs.results)
  #Adjusted P values (FDR Q values)
pdf(file="volcano2.pdf")
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))
with(subset(topT , padj<.05), points(log2FoldChange, -log10(padj), pch=20,col="black"))
with(subset(topT, padj<0.05 & (log2FoldChange)>0), points(log2FoldChange, -log10(padj), pch=20, col="violetred",cex=0.5))
with(subset(topT, padj<0.05 & (log2FoldChange)<(0*-0)), points(log2FoldChange, -log10(padj), pch=20, col="royalblue4",cex=0.5))
dev.off()

#save list of DEGs
DEGs_mangrove_vs_reef <- as.data.frame(subset(DEGs.results, padj<0.05))
DEGs_mangrove_vs_reef_ordered <- order(DEGs_mangrove_vs_reef$padj) #Order p-values by smallest value first
DEGs_mangrove_vs_reef$contrast <- as.factor(c("mangrove_vs_reef"))
DEGs_mangrove_vs_reef$gene_id  <- rownames(DEGs_mangrove_vs_reef)
rownames(DEGs_mangrove_vs_reef) <- NULL
write.csv(DEGs_mangrove_vs_reef, "DEGs_mangrove_vs_reef.csv")
```
#### Plot DEGs 
```{r}
#identify signficant pvalues with 5%FDR
DEGs.results$gene_id  <- rownames(gdds)
sig <- subset(DEGs.results, padj<0.05,)
rownames(sig) <- sig[,7] #rename rownames of sig as column 7
#subset list of sig transcripts from original count data
sig.list <- gdds[which(rownames(gdds) %in% rownames(sig)),]

#apply a vst transformation to minimize effects of small counts and normalize by library size
SF.gdds <- estimateSizeFactors(gdds)
print(sizeFactors(SF.gdds))
vstsig<-vst(sig.list, blind=FALSE)

PCA.sig.degs <- plotPCA(rsig, intgroup="location")+geom_text(aes(label=name),vjust=2) #Plot PCA of all samples for DEG only
PCA.sig.degs #view plot

pdf(file="PCA_sigDEGs_mangrove_vs_reef.pdf")
PCA.sig.degs #view plot
dev.off()

## Heatmap 
#make an expression object
#difference in expression compared to average across all samples
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(MetBrewer)

colors = met.brewer("OKeeffe1",n=9,type="continuous", direction=-1)

Sym.mat <- assay(vstsig)
Sym.mat <- Sym.mat - rowMeans(Sym.mat)
Sym.df <- data.frame(colData(rsig)[c("location")])
pdf(file="SigDEGs_Heatmap_adults_changeColors.pdf")
pheatmap(Sym.mat, annotation_col = Sym.df, clustering_method = "average",
         clustering_distance_rows="euclidean", 
         #color = viridis(400),
         color = colors,
         show_rownames =FALSE, cluster_cols=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

dev.off()
```

#### Visualize differentially-expressed genes: CLUSTERS
```{r}
##### Subset and Log-transform the count data
#Subset the gene count matrix by the list of DEGs
DEGs.results$gene_id  <- rownames(gdds)
sig <- subset(DEGs.results, padj<0.05,)
rownames(sig) <- sig[,7]
sig.list <- gdds[which(rownames(gdds) %in% rownames(sig)),]

#apply a rlog transformation to minimize effects of small counts and normalize wrt library size
vstsig <- vst(sig.list, blind=FALSE)

#Make a matrix for computing similarity
mat <- assay(vstsig) #make an expression object
mat <- mat - rowMeans(mat) #difference in expression compared to average across all samples

### Compute the optimal number of clusters for plotting
#Find the optimum number of clusters using 30 indexes with the NbClust() package. 

nb <- NbClust(mat, distance = "euclidean", min.nc = 2,
              max.nc = 5, method = "kmeans")
#[1]According to the majority rule, the best number of clusters is  2 
 
fviz_nbclust(nb)

calc.kmeans <- kmeans(mat, 2)
cluster_res <- data.frame(gene_id = names(calc.kmeans$cluster), cluster = calc.kmeans$cluster)
# 

DEGs <- as.data.frame(subset(DEGs.results, padj<0.05))
DEGs$contrast <- as.factor(c("mangrove_vs_reef"))

DEGs_cluster <- merge(DEGs, cluster_res, by = "gene_id")
write.csv(DEGs_cluster, "DEGs_cluster.csv")

#### Plot a heatmap of differentially-expressed genes

DEGs_cluster_heat <- read.csv("DEGs_clusters.csv", header = TRUE, sep = ",")[,-c(1)]
DEGs_cluster_heat <- subset(DEGs_cluster, select = c(gene_id, cluster))

#Prepare annotations
hm_ann_row <- unique(DEGs_cluster_heat)
rownames(hm_ann_row) <- hm_ann_row$gene_id
hm_ann_row <- subset(hm_ann_row, select=cluster)
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)
hm_ann_row <- as.matrix(hm_ann_row[rownames(mat),])
hmlocation <- colData(gvst)[c("location")]

hmlocation$location <- factor(hmlocation$location, levels=c("reef", "mangrove"))
hm_ann_col <- HeatmapAnnotation(df=hmlocation, col = list(depth=c("mangrove" ="cadetblue", "reef"  ="indianred3"))) #make dataframe for column naming

library("RColorBrewer")
library(MetBrewer)

DEGs_clusterHeatmap <-  Heatmap(mat, column_title = "Location", 
                                      name = "expression",
                                      #col=brewer.pal(11,"RdBu"),
                                      #col= met.brewer("Hiroshige",type="continuous"),
                                      col= met.brewer("Benedictus",direction=-1),
                                     #scale_fill_gradientn(colors= met.brewer("Hiroshige",type="continuous")),
                                      show_row_names = FALSE, top_annotation = hm_ann_col, show_column_names = FALSE, row_dend_side = "left" ,
                                      column_split = 2, column_dend_height = unit(0.5, "in"),
                                      km = 2, row_km_repeats = 100, row_title = c("Cluster1", "Cluster2"),
                                      row_gap = unit(2.5, "mm"), border = TRUE,
                                      column_names_gp =  gpar(fontsize = 10)); DEGs_clusterHeatmap
                                
png("DEGs_clusterHeatmap.png")
DEGs_clusterHeatmap
dev.off()

pdf(file="DEGs_clusterHeatmap.pdf")
DEGs_clusterHeatmap #view plot
dev.off()
```

