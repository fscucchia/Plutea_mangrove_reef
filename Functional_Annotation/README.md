
**Functional annotation of *Porites lutea* reference genome**

# Create and activate a conda environment
```
conda create -n functional_annotation
conda activate functional_annotation
```

# Install all necessary programs within your conda environment 

- [DIAMOND](http://www.diamondsearch.org/) Search Program
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
#added the path in bashrc
```
### 1) Find homologous sequences

#### i) Download nr database
Download the nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). Then, use Diamond's ```makedb``` command to format the database in a Diamond-friendly format. 

```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz #download nr database in fasta format, I downloaded version of 2021-11-27 
diamond makedb --in nr.gz -d nr     
```
#### ii) Run DIAMOND Search

As input, DIAMOND requires your reference sequences (either protein or CDS nucleotides), and a path to your nr database. 

- Download *P. lutea* reference sequences from reefgenomics
```
wget http://plut.reefgenomics.org/download/plut2v1.1.proteins.fasta
```
- Run sequence alignment against the nr database
Run the script [`diamond_nrNCBI_PluteaCDS.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/edit/main/Functional_Annotation/diamond_nrNCBI_PluteaCDS.sh) that contains the following commands: 
```
diamond blastp -d /data/home/databases/nr_NCBI/nr.dmnd -q /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.proteins.fasta -o $O/Plut.diamondBlastpNCBInr_20012022 -f 6 -b 20 --more-sensitive -e 0.00001 -k1       
#Use -f 6 to get the diamond results in tab format
#Reported 28844 pairwise alignments with diamond blastp, took 14 hours
```
#### iii) Run BLASTp Search
**Blastp: Align _P. lutea_ protein query sequences against a protein reference database (NCBI nr)**
Run the script [`Blastp_nrNCBI_Plutea.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/edit/main/Functional_Annotation/Blastp_nrNCBI_Plutea.sh) that contains the following commands: 

- Prepare the nr database for Blastp
```
gunzip -c nr.gz     #I renamed the output file adding ".fasta"
makeblastdb -in nr.fasta -dbtype prot #took 10 hours
```
- Run blastp on nr database 
```
$ blastp -query $Q/plut2v1.1.proteins.fasta -db $H/nr.fasta -num_threads 20 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out $O/Plut.annot_blastp.200120222
```

### 2) Uniprot 

#### Annotate sequences with GO terms using the UniProt-GOA database and results from blastp and diamond 
With excel, I selected from the results of diamond blastp only the column with the protein accession IDs (second column). Used this list called `Plutea_accession_IDs_DIAMONDblastp.txt` as input into [Uniprot mapping services](https://www.uniprot.org/uploadlists/) to map to the UniProtKB and UniParc IDs databases (contains GO and KO info).

- 381 out of 21664 RefSeq Protein identifiers were successfully mapped to 381 UniProtKB IDs -> Saved the results as txt `Plutea_Uniprot_RefSeqProtein_to_UniprotKB_IDs_diamondBlastp.tab`.

- 8 out of 21664 UniProtKB AC/ID identifiers were successfully mapped to 8 UniProtKB/Swiss-Prot IDs -> Saved the results as txt `Plutea_Uniprot_UniProtKB_to_UniProtKB/SwissProt_diamondBlastp.tab`.

- 1 out of 21664 EMBL/GenBank/DDBJ identifiers were successfully mapped to 1 UniProtKB/Swiss-Prot ID -> Saved the results as txt `Plutea_Uniprot_EMBLGenBankDDBJ_to_UniProtKBSwiss-Prot_diamondBlastp.tab`.

- 2472 out of 21664 EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to 2472 UniProtKB IDs -> Saved the results as txt `Plutea_Uniprot_EMBLGenBankDDBJCDS_to_UniProtKB_diamondBlastp.tab`.

- 4 out of 21664 Ensembl Genomes Protein identifiers were successfully mapped to 4 UniProtKB IDs -> Saved the results as txt `Plutea_Uniprot_EnsemblGenomesProtein_to_UniProtKB_diamondBlastp.tab`.

### 3) EggNog
EggNOG-mapper is a tool for fast functional annotation of novel sequences. It uses precomputed Orthologous Groups (OGs) and phylogenies from the [EggNOG database](http://eggnog5.embl.de) to transfer functional information from fine-grained orthologs only. Common uses of eggNOG-mapper include the annotation of novel genomes, transcriptomes or even metagenomic gene catalogs.

- Used the [eggnog-mapper](http://eggnog-mapper.embl.de/) to annotate results from previous steps.

### 4) Compilation of the output of different methods (Uniprot and EggNog)

Done in RStudio.
### Set up workspace in R Studio
```{r}
# Load libraries
library("tidyverse") 
library("dplyr")
library("readr")
library("tidyr")
```
#### Diamond_Blastp results
```{r}
blast <- read_tsv("Plut.diamondBlastpNCBInr_20012022.txt", col_names = FALSE)
colnames(blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(blast)
dim(blast)
```

#### Unitprot results
```{r}
#381 out of 21664 RefSeq Protein identifiers were successfully mapped to 381 UniProtKB IDs 
#-> Saved the results as txt `Plutea_Uniprot_RefSeqProtein_to_UniprotKB_IDs_diamondBlastp.tab`.
u1 <- read_tsv("Plutea_Uniprot_RefSeqProtein_to_UniprotKB_IDs_diamondBlastp.tab", col_names = TRUE)
u1 <- u1[,c(1:2,4:10)]
colnames(u1) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")
head(u1)
dim(u1)

#8 out of 21664 UniProtKB AC/ID identifiers were successfully mapped to 8 UniProtKB/Swiss-Prot IDs 
#-> Saved the results as txt `Plutea_Uniprot_UniProtKB_to_UniProtKB/SwissProt_diamondBlastp.tab`
u2 <- read_tsv("Plutea_Uniprot_UniProtKB_to_UniProtKBSwissProt_diamondBlastp.tab", col_names = TRUE)
u2 <- u2[,c(1:2,4:10)]
colnames(u2) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")

#1 out of 21664 EMBL/GenBank/DDBJ identifiers were successfully mapped to 1 UniProtKB/Swiss-Prot ID 
#-> Saved the results as txt `Plutea_Uniprot_EMBLGenBankDDBJ_to_UniProtKBSwiss-Prot_diamondBlastp.tab`.
u3 <- read_tsv("Plutea_Uniprot_EMBLGenBankDDBJ_to_UniProtKBSwiss-Prot_diamondBlastp.tab", col_names = TRUE)
u3 <- u3[,c(1:2,4:10)]
colnames(u3) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")

#2472 out of 21664 EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to 2472 UniProtKB IDs 
#-> Saved the results as txt `Plutea_Uniprot_EMBLGenBankDDBJCDS_to_UniProtKB_diamondBlastp.tab`.
u4 <- read_tsv("Plutea_Uniprot_EMBLGenBankDDBJCDS_to_UniProtKB_diamondBlastp.tab", col_names = TRUE)
u4 <- u4[,c(1:2,4:10)]
colnames(u4) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")

#4 out of 21664 Ensembl Genomes Protein identifiers were successfully mapped to 4 UniProtKB IDs 
#-> Saved the results as txt `Plutea_Uniprot_EnsemblGenomesProtein_to_UniProtKB_diamondBlastp.tab`.
u5 <- read_tsv("Plutea_Uniprot_EnsemblGenomesProtein_to_UniProtKB_diamondBlastp.tab", col_names = TRUE)
u5 <- u5[,c(1:2,4:10)]
colnames(u5) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")

#15150 out of 21760 RefSeq Protein identifiers were successfully mapped to 15150 UniParc IDs
#-> Saved the results as txt `Plutea_Uniprot_RefSeqProtein_to_UniParcIDs.tab` (this is from the blastx of 
#diamond, the apping RefSeq to Uniparc doesn't work anymore). Then I mapped the Uniparc to UniProtKB to get 
#the Gos. Mapping to Uniparc doesn't give GOs.
u6a <- read_tsv("Plutea_Uniprot_RefSeqProtein_to_UniParcIDs.tab", col_names = TRUE)
u6b <- read_tsv("GO_Uniparc_to_Uniprot.tab", col_names = TRUE)

u6 <- merge(u6a, u6b, by = "UPA", all = TRUE)
u6 <- na.omit(u6, cols = c(8:16)) 
#write.csv(u6, file = "RefSeq_Uniparc_withGO.csv", row.names = FALSE)
u6 <- u6[,c(2,8,10:16)]
colnames(u6) <- c("top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids")
```
#### Compile Uniprot results
```{r}
Uniprot_results <- bind_rows(u1, u2, u3, u4, u5, u6)
Uniprot_results <- unique(Uniprot_results)
Uniprot_results$go_ids <- gsub(" ", "", Uniprot_results$go_ids)

nrow(filter(Uniprot_results, grepl("GO",go_ids))) #Genes with GO terms
```

#### EggNog results
```{r}
egg <- read_tsv("EggNog_adjusted.txt", col_names = TRUE) #the top_hit names re problematic, so I used the suffix way below to separate the names 

suffix1 <- read_tsv("suffix1.txt", col_names = TRUE)
suffix1$seed_ortholog <- sub("$", ".1", suffix1$seed_ortholog)

suffix2 <- read_tsv("suffix2.txt", col_names = TRUE)
suffix2$seed_ortholog <- sub("$", ".2", suffix2$seed_ortholog)

suffix_final <- bind_rows(suffix1, suffix2)
suffix_final <- suffix_final[,-3] 

notsuffix <- read_tsv("notsuffix.txt", col_names = TRUE)

egg_final <- bind_rows(suffix_final, notsuffix)

colnames(egg_final) <- c("seqName", "top_hit", "evalue", "score", "OGs", "max_annot", "COG_cat", "description", "go_ids")
nrow(filter(egg_final, grepl("GO",go_ids))) #Genes with GO terms... Commented out because all have go terms
```
#### Find unique and overlapping GO terms
```{r}
#Generate lists of GO terms for each method
Uniprot_GO <- select(Uniprot_results, top_hit, go_ids)
Uniprot_GO_splitted <- strsplit(as.character(Uniprot_GO$go_ids), ";") #split into multiple GO ids
gene_ontology <- data.frame(v1 = rep.int(Uniprot_GO$top_hit, sapply(Uniprot_GO_splitted, length)), v2 = unlist(Uniprot_GO_splitted)) #list all genes with each of their GO terms in a single row
colnames(gene_ontology) <- c("prot_id", "GO.ID")

gene_ontology <- separate(gene_ontology, GO.ID, into = c("GO.ID", "ontology", "term"), sep=" ") #Split GO.ID, terms and ontologies into separate columns
Uniprot.GOterms <- select(gene_ontology, prot_id, GO.ID)

Uniprot.GOterms$GO.ID<- as.character(Uniprot.GOterms$GO.ID)
Uniprot.GOterms[Uniprot.GOterms == 0] <- "unknown"
Uniprot.GOterms$GO.ID <- replace_na(Uniprot.GOterms$GO.ID, "unknown")
Uniprot.GOterms$GO.ID <- as.factor(Uniprot.GOterms$GO.ID)
Uniprot.GOterms$prot_id <- as.factor(Uniprot.GOterms$prot_id)
Uniprot.GOterms$GO.ID <- gsub(" ", "", Uniprot.GOterms$GO.ID)
Uniprot.GOterms <- unique(Uniprot.GOterms)
nrow(Uniprot.GOterms)

egg_GO <- select(egg_final, top_hit, go_ids)
egg_splitted <- strsplit(as.character(egg_GO$go_ids), ";") #split into multiple GO ids
egg_gene_ontology <- data.frame(v1 = rep.int(egg_GO$top_hit, sapply(egg_splitted, length)), v2 = unlist(egg_splitted)) #list all genes with each of their GO terms in a single row
colnames(egg_gene_ontology) <- c("prot_id", "GO.ID")
egg_gene_ontology <- separate(egg_gene_ontology, GO.ID, into = c("GO.ID", "ontology", "term"), sep=" ") #Split GO.ID, terms and ontologies into separate columns
egg_GOterms <- select(egg_gene_ontology, prot_id, GO.ID)
egg_GOterms$GO.ID<- as.character(egg_GOterms$GO.ID)
egg_GOterms[egg_GOterms == 0] <- "unknown"
egg_GOterms$GO.ID <- replace_na(egg_GOterms$GO.ID, "unknown")
egg_GOterms[egg_GOterms == "-"] <- "unknown"
egg_GOterms$GO.ID <- as.factor(egg_GOterms$GO.ID)
egg_GOterms$prot_id <- as.factor(egg_GOterms$prot_id)
egg_GOterms$GO.ID <- gsub(" ", "", egg_GOterms$GO.ID)
egg_GOterms <- unique(egg_GOterms)
nrow(egg_GOterms)
```

#### Find intersections and unique results for each methods
```{r}
UE <- intersect(egg_GOterms, Uniprot.GOterms) #EggNog and Uniprot intersection
nrow(UE)
Uunique <- setdiff(Uniprot.GOterms, egg_GOterms) #Uniprot unique
nrow(Uunique)
Eunique <- setdiff(egg_GOterms, Uniprot.GOterms) #EggNog unique
nrow(Eunique)
```

#### Merge Annotations
```{r}
# Match top_hits with description
Plutea_annot <- left_join(blast, egg_final, by="seqName")
Plutea_annot <- select(Plutea_annot, seqName, top_hit.x, length, evalue.x, bitscore, go_ids)
Plutea_annot <- rename(Plutea_annot, "top_hit"="top_hit.x")
Plutea_annot <- left_join(Plutea_annot, Uniprot_results, by="top_hit")
Plutea_annot$GO <- paste(Plutea_annot$go_ids.x, Plutea_annot$go_ids.y, sep=';') #generate new column with concatenated GO IDs
#Plutea_annot$GO_terms <- paste(Plutea_annot$GO_names, Mcap_annot$gene_ontology, sep=';') #generate new column with concatenated GO IDs
Plutea_annot2 <- select(Plutea_annot,-c("go_ids.x", "go_ids.y", "length.y"))
colnames(Plutea_annot2) <- c("prot_id", "top_hit", "length","eValue", "bitscore","UniProtKB_entry", "status", "protein_names", "gene_names","organism","GO_terms","GO_IDs")

write_tsv(Plutea_annot2, "Plutea_annot_GO.tsv")
```


