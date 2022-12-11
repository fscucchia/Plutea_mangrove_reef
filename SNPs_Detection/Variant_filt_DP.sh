#!/bin/sh

F="/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads"
REF="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"
OUT2="GVCFall_morefilter"
O="/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration"
###
#DP_MIN=2
DP_MIN=20.000
###

gatk VariantFiltration --reference "${REF}" --variant "${OUT2}_SNPs_VarScores_filterPASSED.vcf" --output "${OUT2}_SNPs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT2}_SNPs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1
gatk VariantFiltration --reference "${REF}" --variant "${OUT2}_INDELs_VarScores_filterPASSED.vcf" --output "${OUT2}_INDELs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT2}_INDELs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1