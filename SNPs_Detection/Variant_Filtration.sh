
#!/bin/sh

F="/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads"
REF="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"
OUT2="GVCFall_morefilter"
O="/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration"

####
SNP_QD_MIN=20.000
SNP_MQ_MIN=50.000
SNP_FS_MAX=5.000
#SNP_SOR_MAX=2.500
SNP_SOR_MAX=10.000

INDEL_QD_MIN=20.000
INDEL_MQ_MIN=45.000
INDEL_FS_MAX=5.000
INDEL_SOR_MAX=10.000

#DP_MIN=2
DP_MIN=20.000
####

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$F --job-name "SelectVariants" -o "$O/SelectVariants.out" -e "$O/SelectVariants.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk SelectVariants --reference "${REF}" --variant "cohort_genotypes.vcf.gz" --output "/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration/${OUT2}_SNPs.vcf.gz"   -select-type SNP 1> "/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration/${OUT2}_SNPs.vcf.gz.log" 2>&1
             gatk SelectVariants --reference "${REF}" --variant "cohort_genotypes.vcf.gz" --output "/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration/${OUT2}_INDELs.vcf.gz" -select-type INDEL 1> "/data/home/Plutea_mangroves/output_Plutea_host/Variant_Filtration/${OUT2}_INDELs.vcf.gz.log" 2>&1
            "

elif [ $1 -eq 2 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_SNPs.vcf.gz"   --output "${OUT2}_SNPs.table"   -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum 1> "${OUT2}_SNPs.table.log" 2>&1
              gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_INDELs.vcf.gz" --output "${OUT2}_INDELs.table" -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum 1> "${OUT2}_INDELs.table.log" 2>&1
              "

elif [ $1 -eq 3 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "variants_scores" -o "$O/variants_scores.out" -e "$O/variants_scores.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap ". $CONDA; conda activate R.settings2
                        Rscript plot_variants_scores.R "${OUT2}.variants_scores" "${OUT2}_SNPs.table" "${OUT2}_INDELs.table" \
	                   $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
	                   $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
	                   1> "${OUT2}.variants_scores.log" 2>&1
                        conda deactivate
                        "

elif [ $1 -eq 4 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantFiltration" -o "$O/VariantFiltration.out" -e "$O/VariantFiltration.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "chmod u+x Variant_filt.sh; ./Variant_filt.sh"

elif [ $1 -eq 5 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_SNPs_VarScores_filterPASSED.vcf"   --output "${OUT2}_SNPs_VarScores_filterPASSED.table" \
	         -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
	         1> "${OUT2}_SNPs_VarScores_filterPASSED.table.log" 2>&1
              gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_INDELs_VarScores_filterPASSED.vcf" --output "${OUT2}_INDELs_VarScores_filterPASSED.table" \
	         -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
	         1> "${OUT2}_INDELs_VarScores_filterPASSED.table.log" 2>&1
              . $CONDA; conda activate R.settings2
                        Rscript check_filtering.R "${OUT2}_SNPs_VarScores_filterPASSED.table" "${OUT2}_INDELs_VarScores_filterPASSED.table" \
                        $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
                        $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
                    	1> "${OUT2}.check_filtering.log" 2>&1
                        Rscript plot_variants_scores2.R "${OUT2}.variants_scores_afterFiltering" "${OUT2}_SNPs_VarScores_filterPASSED.table" "${OUT2}_INDELs_VarScores_filterPASSED.table" \
                        $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
                        $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
	                   1> "${OUT2}.variants_scores_afterFiltering.log" 2>&1
                        conda deactivate"

elif [ $1 -eq 6 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads/cohort_genotypes.vcf.gz" --output "${OUT2}.DP.table" -F CHROM -F POS -GF GT -GF DP 1> "${OUT2}.DP.table.log" 2>&1
              "

# elif [ $1 -eq 7 ]; then
#      sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "Extract_DP" -o "$O/Extract_DP.out" -e "$O/Extract_DP.err" \
#      -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
#   	 --wrap "chmod u+x Extract_DP.sh; ./Extract_DP.sh"
              

elif [ $1 -eq 7 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "Variant_filt_DP" -o "$O/Variant_filt_DP.out" -e "$O/Variant_filt_DP.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "chmod u+x Variant_filt_DP.sh; ./Variant_filt_DP.sh"

elif [ $1 -eq 8 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk SelectVariants --reference "${REF}" --variant "${OUT2}_SNPs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1
     gatk SelectVariants --reference "${REF}" --variant "${OUT2}_INDELs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT2}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT2}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1"

elif [ $1 -eq 9 ]; then
     sbatch -N1 --workdir=$O --job-name "Extract_plot_DP" -o "$O/Extract_plot_DP.out" -e "$O/Extract_plot_DP.err" \
     -p hiveunlim \
  	 --wrap "chmod u+x Extract_plot_DP_info.sh; ./Extract_plot_DP_info.sh"

elif [ $1 -eq 10 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	--wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" -F CHROM -F POS -GF GT -GF AD -GF DP \
	1> "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table.log" 2>&1"

elif [ $1 -eq 11 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf"   --output "${OUT2}_SNPs_filterPASSED_final.table"   -F CHROM -F POS -GF GT \
	1> "${OUT2}_SNPs_filterPASSED_final.table.log" 2>&1
     gatk VariantsToTable --reference "${REF}" --variant "${OUT2}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT2}_INDELs_filterPASSED_final.table" -F CHROM -F POS -GF GT \
	1> "${OUT2}_INDELs_filterPASSED_final.table.log" 2>&1"

     
fi


