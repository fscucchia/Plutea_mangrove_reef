H="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
F="/data/home/Plutea_mangroves/scripts"
G="/data/home/databases/Plutea_genome_reefgenomics"

########################

if [ $1 -eq 1 ]; then
   mkdir -p $H
   sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$H --job-name "symblink" -o "$H/symblink.out" -e "$H/symblink.err" \
   -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
	 --wrap ". $CONDA; conda activate newrnapipeline; ln -s /data/home/Plutea_mangroves/output_Plutea_host/fastqc_filtered/*gz.filtered ./; conda deactivate"

elif [ $1 -eq 2 ]; then
     mkdir -p $G
     sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$G --job-name "hisat" -o "$H/hisat.out" -e "$H/hisat.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap ". $CONDA; conda activate newrnapipeline; hisat2-build -f /data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta ./Plut_ref
	                   mkdir -p $F
					   chmod u+x StyHISAT_withSummary.sh
					   ./StyHISAT_withSummary.sh 
					   conda deactivate"

elif [ $1 -eq 3 ]; then # multiqc
		mkdir -p $H
		sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$H --job-name "multiqc" -o "$H/multiqc_hisat".out -e "$H/multiqc_hisat".err \
		-p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
		--wrap ". $CONDA; conda activate multiqc_env; mkdir -p $H/multiqc_hisat; multiqc $H; conda deactivate"		

fi	
						