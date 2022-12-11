

###################### Genetic differentiation analisys using SNPRelate and HIERFSTAT 

setwd("/data/home/Plutea_mangroves/output_Plutea_host/genetic_diff")

## Load libraries

library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)
library(ggplot2)

## Load the pruned SNPs (pruned with plink)
bed.fn <- "/data/home/Plutea_mangroves/output_Plutea_host/genotype_analysis/P_lutea_vcf_pruned.bed"
fam.fn <- "/data/home/Plutea_mangroves/output_Plutea_host/genotype_analysis/P_lutea_vcf_pruned.fam"
bim.fn <- "/data/home/Plutea_mangroves/output_Plutea_host/genotype_analysis/P_lutea_vcf_pruned.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned.gds")
snpgdsSummary("plink_pruned.gds")
# [1]
# The file name: /data/home/Plutea_mangroves/output_Plutea_host/SNPRelate/plink_pruned.gds
# The total number of samples: 8
# The total number of SNPs: 22925
# SNP genotypes are stored in SNP-major mode (Sample X SNP).

## Open the gds file
genofile <- snpgdsOpen("plink_pruned.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.id

## Get population information
pop_code <- scan("/data/home/Plutea_mangroves/output_Plutea_host/SNPRelate/pop.txt", what=character())
#see the file here https://github.com/fscucchia/Plutea_mangrove_reef/blob/main/SNPs_Detection/pop.txt

head(cbind(sample.id, pop_code)) #assumes the order of sample IDs is as the same as population codes

#### Identity-By-State Analysis

ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE)
# Identity-By-State (IBS) analysis on genotypes:
# Excluding 110 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# Working space: 8 samples, 26,201 SNPs
#     using 2 (CPU) cores
# IBS:    the sum of all selected genotypes (0,1,2) = 80725

#To perform multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances:
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

pdf(paste0('MDS_IBS','.pdf'))
plot(x, y, col=race, xlab = "", ylab = "",
    main = "Multidimensional Scaling Analysis (IBS)2")
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))
dev.off()

#To perform cluster analysis on the n×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score:
set.seed(999)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE))

## Draw philogenetic tree
snpHCluster =  snpgdsHCluster((snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE)),sample.id=NULL, need.mat=TRUE, hang=0.01)
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)
pdf(paste0('Tree','.pdf'), width=15, height=20)
snpgdsDrawTree(cutTree, main = "Phylogenetic Tree",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular")
dev.off()

######### Determine groups of individuals automatically https://rdrr.io/bioc/SNPRelate/man/snpgdsCutTree.html
rv <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

rv1 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code),z.threshold=15, n.perm = 5000, 
    col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=4,
    label.H=FALSE, label.Z=TRUE, verbose=TRUE)


# Determine groups by permutation (Z threshold: 15, outlier threshold: 5):
# Create 1 groups.
pdf(paste0('dendrogram1','.pdf'))
plot(rv1$dendrogram, leaflab="none", main="HapMap Phase II")
dev.off()

############# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

############# cluster individuals by Z score, specifying 'clust.count'  
pdf(paste0('dendrogram3','.pdf'))
snpgdsDrawTree(rv2, rv2$clust.count, main="HapMap Phase II",
    edgePar = list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),
    labels = c("mangrove", "reef"), y.label=0.1,
    y.label.kinship=TRUE)
dev.off()


################# Compute Fst statistic

x2<-2-snpgdsGetGeno(genofile) 
#[1]Genotype matrix: 8 samples X 22925 SNPs    

#add column with samples names as numbers
vec <- c("R16b","TM7","T_M2","T_M4C","T_M5","T_R11B",
          "T_R15B","T_R9D") 

locus_transpose1 <- cbind(x2, sample = vec)  

locus_transpose2=as.data.frame(locus_transpose1)
colnames(locus_transpose2[1:10])
rownames(locus_transpose2[1:10])

#move last column with samples names to first
locus_transpose3 <- locus_transpose2 %>%
  select(sample, everything())
head(locus_transpose3[,1:10])

#add column with with groups as numbers
vec2 <- c("1","2","2","2","2","1","1","1")

locus_transpose4 <- cbind(locus_transpose3, group = vec2) 

#move last columns with depth and colony to 2nd and 3rd places
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#move group column after sample column using the moveme function above
locus_transpose5=locus_transpose4[moveme(names(locus_transpose4), "group after sample")]
head(locus_transpose5[,1:10])

locus_transpose5$sample <- NULL #don't need this column

## Fst with package hierfstat, Method: Weir & Cockerham, 1984

library(hierfstat)

locus_transpose6 <- lapply(locus_transpose5,as.numeric) #pairwise.WCfst only works with numerics
locus_transpose7=as.data.frame(locus_transpose6)

Fst_between_groups = pairwise.WCfst(locus_transpose7,diploid=TRUE) #Fst between mangrove and reef sites
Fst_between_groups
#[1] Fst_between_groups  
#            1          2
# 1         NA 0.05110393
# 2 0.05110393         NA

