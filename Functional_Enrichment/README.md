
## GO enrichment analysis _P. lutea_ adult corals from mangrove and reed sites ######

```
#This script is based on the work of Erin Chille (https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/tree/main/5-Planula-GO-Enrichment-Analysis) with some modifications and adjustments.

### Set up workspace

rm(list = ls()) #clear environment

#Load libraries

library(goseq)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(patchwork)
library(dplyr)
library(tidyr)

#treatment information
treatmentinfo <- read.csv("RNAseq_data.csv", header = TRUE, sep = ";")
str(treatmentinfo)
head(treatmentinfo)

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id"))
gcount$gene_id <- rownames(gcount)

### Load DESeq2 results
DEG.res <- read.csv("DEGs_clusters.csv")[,-1]
nrow(DEG.res)

### Load annotations
annot <- read_tsv( "Plutea_annot_GO_cleaned.tsv.txt", col_names = TRUE) #biological annotation information
#remove duplicate genes
annot = annot[order(annot[,'prot_id'],-annot[,'length']),]
annot = annot[!duplicated(annot$prot_id),]

Go.ref <- subset(annot, select= c(prot_id, length)) #Select only relevant information

names(Go.ref)[1] <- "gene_id" #rename column
Go.ref <- unique(Go.ref)

#Filter gcount by available annotations
Go.ref_gcount_merged <- merge(gcount, Go.ref, by = "gene_id")

#Make a dataframe containing the gene_ids and cluster for each cluster.
#Select only gene_id and cluster from DEseq2 res
DEGclust <- subset(DEG.res, select=c(gene_id, cluster))
DEGclust <- unique(DEGclust)
clust1 <- filter(DEGclust, cluster=="1")
nrow(clust1) #nrow clust1
clust2 <- filter(DEGclust, cluster=="2")
nrow(clust2) #nrow clust2

#Set ID and gene length vectors, and make a binary matrix indicating which genes are differentially expressed. These are used as input to nullp, which for calculates a Probability Weighting Function for each set of DEGs.
#Make ID and length vectors
Go.ref_gcount_merged <- unique(Go.ref_gcount_merged)
dim(Go.ref_gcount_merged)
IDvector <- Go.ref_gcount_merged$gene_id
lengthVector <- Go.ref_gcount_merged$length

#Cluster 1
clust1.genes <- as.vector(clust1$gene_id)
clust1.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust1.genes)
names(clust1.genes)=Go.ref_gcount_merged$gene_id
length(clust1.genes)
length(names(clust1.genes))
length(unique(names(clust1.genes)))

#Cluster 2
clust2.genes <- as.vector(clust2$gene_id)
clust2.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust2.genes)
names(clust2.genes)=Go.ref_gcount_merged$gene_id
length(clust2.genes)
length(names(clust2.genes))
length(unique(names(clust2.genes)))

pwf.C1<-nullp(DEgenes=clust1.genes, id=IDvector, bias.data=lengthVector) #weight vector by length of gene
pwf.C2<-nullp(clust2.genes, IDvector, bias.data=lengthVector) #weight vector by length of gene

#### Prepare GO term dataframe

GO.annot <- subset(annot, select= c(prot_id, GO_IDs)) #Select only relevant information

names(GO.annot)[1] <- "gene_id" #rename column

GO.annot.na <- filter(GO.annot, GO_IDs!="NA;NA") #Remove NAs
GO.annot.na_cleaned <- GO.annot.na[!grepl("-", GO.annot.na$GO_IDs), ] #Remove "-", that is genes with no GO term

splitted <- strsplit(as.character(GO.annot.na_cleaned$GO_IDs), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot.na_cleaned$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")

GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
GO.terms <- unique(GO.terms)

nrow(GO.terms)/length(unique(GO.terms$gene_id))
#[1] 88.04816

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.C1 <- goseq(pwf.C1, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GOwall.C2 <- goseq(pwf.C2, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#Find only enriched GO terms that are statistically significant at cutoff
C1.GO.sigp<-GOwall.C1$category[GOwall.C1$over_represented_pvalue<.05]
C1.GO.sigp<-data.frame(C1.GO.sigp)
colnames(C1.GO.sigp) <- c("category")
C1.GO.sigp <- merge(C1.GO.sigp, GOwall.C1, by="category")
C1.GO.sigp <- C1.GO.sigp[order(C1.GO.sigp$ontology, C1.GO.sigp$over_represented_pvalue, -C1.GO.sigp$numDEInCat),]
C1.GO.sigp$term <- as.factor(C1.GO.sigp$term)

nrow(C1.GO.sigp) 

C2.GO.sigp<-GOwall.C2$category[GOwall.C2$over_represented_pvalue<.05]
C2.GO.sigp<-data.frame(C2.GO.sigp)
colnames(C2.GO.sigp) <- c("category")
C2.GO.sigp <- merge(C2.GO.sigp, GOwall.C2, by="category")
C2.GO.sigp <- C2.GO.sigp[order(C2.GO.sigp$ontology, C2.GO.sigp$over_represented_pvalue, -C2.GO.sigp$numDEInCat),]
C2.GO.sigp$term <- as.factor(C2.GO.sigp$term)

nrow(C2.GO.sigp)

###Save significant terms
write.csv(C1.GO.sigp, "C1.GO.sigp.csv", row.names = FALSE)
write.csv(C2.GO.sigp, file = "C2.GO.sigp.csv", row.names = FALSE)
```

### Find GOslim terms
```

C1.GO.sigp$dir <- "C1 (Down)" #down in reef
C2.GO.sigp$dir <- "C2 (Up)"   #up in reef
all_GO <- bind_rows(C1.GO.sigp, C2.GO.sigp)  #bind rows

#Run GOslim to get broader categories
slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database

## BP
BP_GO <- all_GO %>%
  filter(ontology=="BP")
BPGO_collection <- GOCollection(BP_GO$category) #Make library of query terms
slims_bp <- data.frame(goSlim(BPGO_collection, slim, "BP")) #Find common parent terms to slim down our list
slims_bp$category <- row.names(slims_bp) #save rownames as category

## MF
MF_GO <- all_GO %>%
  filter(ontology=="MF")
MFGO_collection <- GOCollection(MF_GO$category) #Make library of query terms
slims_mf <- data.frame(goSlim(MFGO_collection, slim, "MF")) #Find common parent terms to slim down our list
slims_mf$category <- row.names(slims_mf) #save rownames as category

## CC
CC_GO <- all_GO %>%
  filter(ontology=="CC")
CCGO_collection <- GOCollection(CC_GO$category) #Make library of query terms
slims_cc <- data.frame(goSlim(CCGO_collection, slim, "CC")) #Find common parent terms to slim down our list
slims_cc$category <- row.names(slims_cc) #save rownames as category
#

#Get mapped terms, using functions from Sam White's Biostars [post](https://support.bioconductor.org/p/128407/#128409).
#Write function mappedIds to get the query terms that mapped to the slim categories
mappedIds <-
  function(df, collection, OFFSPRING) #the command to run requires a dataframe of slim terms, like slims_MF above, your list of query terms, and the offspring from the GOCollection by goSlim
  {
    map <- as.list(OFFSPRING)[rownames(df)] # Subset GOcollection offspring by the rownames of your dataframe
    mapped <- lapply(map, intersect, ids(collection)) #Find the terms that intersect between the subset made above of your query terms and the GOids from the GO collection
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L)) #Add column "go_terms" with matching terms 
    df #show resulting dataframe
  }

#Run function for MF and BP terms
BPslim <- mappedIds(slims_bp, BPGO_collection, GOBPOFFSPRING)
MFslim <- mappedIds(slims_mf, MFGO_collection, GOMFOFFSPRING)
CCslim <- mappedIds(slims_cc, CCGO_collection, GOCCOFFSPRING)

#Remove duplicate matches, keeping the broader umbrella term
#BP
BPslim <- filter(BPslim, Count>0 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(BPslim$go_terms), ";") #split into multiple GO ids
BPslimX <- data.frame(Term = rep.int(BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
BPslimX <- merge(BPslimX, BPslim[,c(1,3:4)], by="Term") #Add back counts, term, and category info
BPslimX <- unique(setDT(BPslimX)[order(go_term, -Count)], by = "go_term") #remove duplicate offspring terms, keeping only those that appear in the larger umbrella term (larger Count number)
BPslim <- data.frame(slim_term=BPslimX$Term, slim_cat=BPslimX$category, category=BPslimX$go_term) #rename columns
head(BPslim)

#MF
MFslim <- filter(MFslim, Count>0 & Term!="molecular_function") #filter out empty slims and term "molecular function"
MFsplitted <- strsplit(as.character(MFslim$go_terms), ";") #split into multiple GO ids
MFslimX <- data.frame(Term = rep.int(MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
MFslimX <- merge(MFslimX, MFslim[,c(1,3:4)], by="Term")  #Add back counts, term, and category info
MFslimX <- unique(setDT(MFslimX)[order(go_term, -Count)], by = "go_term")  #remove duplicate offspring terms, keeping only
MFslim <- data.frame(slim_term=MFslimX$Term, slim_cat=MFslimX$category, category=MFslimX$go_term) #rename columns
head(MFslim)

#CC
CCslim <- filter(CCslim, Count>0 & Term!="cellular component") #filter out empty slims and term "molecular function"
CCsplitted <- strsplit(as.character(CCslim$go_terms), ";") #split into multiple GO ids
CCslimX <- data.frame(Term = rep.int(CCslim$Term, sapply(CCsplitted, length)), go_term = unlist(CCsplitted)) #list all
CCslimX <- merge(CCslimX, CCslim[,c(1,3:4)], by="Term")  #Add back counts, term, and category info
CCslimX <- unique(setDT(CCslimX)[order(go_term, -Count)], by = "go_term")  #remove duplicate offspring terms, keeping only
CCslim <- data.frame(slim_term=CCslimX$Term, slim_cat=CCslimX$category, category=CCslimX$go_term) #rename columns
head(CCslim)

#Save slim info with GO enrichment info for heatmap dataframes.
GO.BP <- right_join(BPslim, filter(all_GO, ontology=="BP"), by="category") #add back GO enrichment info for each offspring term
GO.MF <- right_join(MFslim, filter(all_GO, ontology=="MF"), by="category") #add back GO enrichment info for each offspring term
GO.CC <- right_join(CCslim, filter(all_GO, ontology=="CC"), by="category")

GO.BP <- na.omit(GO.BP) 
GO.MF <- na.omit(GO.MF) 

write.csv(GO.BP, "GOs per slim_list_BP.csv")
write.csv(GO.MF, "GOs per slim_list_MF.csv")

######## counting occurrences in data.frame - how many terms per each slimGO and per cluster
occurences_GO.BP_clust1<- filter(GO.BP, dir=="C1 (Down)")
occurences_GO.BP_clust1<-table(unlist(occurences_GO.BP_clust1$slim_term))
write.csv(occurences_GO.BP_clust1, "occurences_GO.BP_clust1.csv")

occurences_GO.BP_clust2<- filter(GO.BP, dir=="C2 (Up)")
occurences_GO.BP_clust2<-table(unlist(occurences_GO.BP_clust2$slim_term))
write.csv(occurences_GO.BP_clust2, "occurences_GO.BP_clust2.csv")

occurences_GO.MF_clust1<- filter(GO.MF, dir=="C1 (Down)")
occurences_GO.MF_clust1<-table(unlist(occurences_GO.MF_clust1$slim_term))
write.csv(occurences_GO.MF_clust1, "occurences_GO.MF_clust1.csv")

occurences_GO.MF_clust2<- filter(GO.MF, dir=="C2 (Up)")
occurences_GO.MF_clust2<-table(unlist(occurences_GO.MF_clust2$slim_term))
write.csv(occurences_GO.MF_clust2, "occurences_GO.MF_clust2.csv")

## Make heatmap
library(org.Hs.eg.db) 
library(MetBrewer)

BPplot <- GO.BP %>% group_by(slim_cat) %>% filter(n()>2) %>% mutate(term = fct_reorder(term, -over_represented_pvalue))%>% filter(numInCat>2) %>% ggplot(aes(x = dir, y = term)) + 
  geom_point(aes(color=over_represented_pvalue, size = 3)) + 
  scale_color_gradientn(colors=met.brewer("VanGogh3",type="continuous"))+
  facet_grid(slim_term ~ ontology, scales = "free_y", labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 11, face = "bold"),
                     strip.text.x = element_text(size = 12, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=15),
                     axis.text = element_text(size = 12), legend.position = "None",
                     plot.margin = unit(c(0,1,0,0.25), "cm")) + theme(strip.background =element_rect(fill="aliceblue"))


MFplot <- GO.MF %>% group_by(slim_cat) %>% filter(n()>2) %>% mutate(term = fct_reorder(term, -over_represented_pvalue)) %>% filter(numInCat>2)%>% ggplot(aes(x = dir, y = term)) + 
  geom_point(aes(color=over_represented_pvalue, size = 3)) + 
  scale_color_gradientn(colors=met.brewer("VanGogh3",type="continuous"))+
  scale_y_discrete(position = "right") +
  facet_grid(slim_term ~ ontology, scales = "free_y", labeller = label_wrap_gen(width = 10, multi_line = TRUE), 
             switch="y") + #Put the y facet strips on the left
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y.left = element_text(angle=0, size = 11, face = "bold"),
                     strip.text.x = element_text(size = 12, face = "bold"),
                     axis.title = element_blank(),
                     axis.text = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 11))+
  theme(strip.background =element_rect(fill="aliceblue"))

enrich_heat <- BPplot + MFplot
enrich_heat <- ggdraw(plot = enrich_heat) + draw_plot_label(c("a)", "b)"), c(0, 0.33), c(1, 1), size = 15)
ggsave("enrich_heat_GOslim.pdf", enrich_heat, width = 30, height = 20, units = c("in"))
```










