# Extreme of today may not save corals from the extreme of tomorrow - The role and risks of coral selective adaptation in the face of climate change 

This electronic notebook provides the scripts employed to analyze gene expression dynamics of the coral _Porites lutea_ inhabiting mangrove and reef sites 
in Woody Isles, Australia.

<p align="center">
<img src="https://github.com/fscucchia/Plutea_mangrove_reef/blob/main/media/Site_map_reduced.png?raw=true" width="700"/>
</p>

### RNA-Seq reads quality filtering and mapping

**[Quality filtering and mapping](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping)** - RNA-Seq reads processing included adapter trimming using Cutadapt v2.6 ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200)) and quality filtering using Trimmomatic v0.39 ([Bolger et al., 2014](https://doi.org/10.1093/bioinformatics/btu170)). Reads were aligned to the [_P. lutea_ genome assembly](http://refuge2020.reefgenomics.org/) using HISAT2 v2.2.1 ([Kim et al., 2019](https://www.nature.com/articles/s41587-019-0201-4)). Transcripts assembly and quantification were performed using Stringtie v2.2.5 ([Pertea et al., 2015](https://www.nature.com/articles/nbt.3122)).

### Coral and algal symbiont species identification

**[Symbiont species identification](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Species_Identification)** - High quality reads were BLASTed using Diamond v2.0.11 ([Buchfink et al., 2021](https://www.nature.com/articles/s41592-021-01101-x)) against the [NCBI](https://www.ncbi.nlm.nih.gov/), [Reefgenomics](http://reefgenomics.org/), [Marinegenomics](https://marinegenomics.oist.jp/gallery) and [UQ eSpace](https://espace.library.uq.edu.au/view/UQ:f1b3a11) proteomes databases of Symbiodiniaceae and stony coral species.

### SNPs characterization

**[Coral host SNPs characterization](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/SNPs_Detection)** - Single nucleotide polymorphisms (SNPs) analysis was conducted using the Genome Analysis Toolkit framework (GATK, v4.2.0; ([McKenna et al., 2010](https://doi.org/10.1101/gr.107524.110))) following the recommended RNA-Seq SNPs practice of the Broad Institute ([(Auwera et al. 2013)](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1110s43)), with necessary adjustments for genotype calling in non-model organisms where variants sites are not known beforehand. HISAT-aligned reads were sorted and marked for duplicates, variant calling was performed with the GATK HaplotypeCaller tool ([McKenna et al., 2010](https://doi.org/10.1101/gr.107524.110)) and genotypes were then jointly called using the GATK GenotypeGVCFs tool. The GATK SelectVariants and VariantFiltration tools were used to filter the joined variant-calling matrices for quality by depth. Filtering for linkage disequilibrium was carried out using PLINK (v2.0, ([Purcell et al. 2007](https://www.cell.com/ajhg/fulltext/S0002-9297(07)61352-4)).
To assess genetic differentiation among mangroves and reef corals, the dissimilarity coefficient and the proportions of shared ancestral alleles were estimated by computing pairwise identity-by-state distances using the SNPRelate R package v1.20.1,69 ([Zheng et al., 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3519454/)). 
In addition, the fixation index (_Fst_)([Weir & Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x)) was estimated using the R package [HIERFSTAT](https://cran.r-project.org/web/packages/hierfstat/index.html) v0.5.10. 

### Differential expression

**[Coral host differential expression](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Differential_Gene_Expression)** - DE analysis was conducted using Bioconductor DEseq2 v1.26.0 ([Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)) in the R environment (v3.6.3) considering a single factor (location).

### Functional annotation

**[Coral host functional annotation](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Functional_Annotation)** - Functional annotation of _P. lutea_ protein sequences was performed using Diamond v2.0.11 ([Buchfink et al., 2021](https://www.nature.com/articles/s41592-021-01101-x)), [Uniprot](https://www.uniprot.org/id-mapping/) and [eggNOG](http://eggnog-mapper.embl.de/).

### Gene ontology enrichment

**[Coral host gene ontology enrichment analysis]https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Functional_Enrichment)** - GO enrichment analysis was performed using the package Goseq (v1.42.0; [Young et al. 2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14)) in the R environment. For the resulting enriched GO terms, slim categories were obtained using the goSlim function of the R package GSEABase (v1.52.1) using the [GOslim generic obo](http://current.geneontology.org/ontology/subsets/goslim_generic.obo) v1.2,65 as reference database .
