H="/data/home/Plutea_mangroves/output_Plutea_host/StringTie"
F="/data/home/Plutea_mangroves/scripts"

########################

if [ $1 -eq 1 ]; then
   sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$H \
   -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
	--wrap ". $CONDA; conda activate newrnapipeline
                     stringtie --merge -p 8 -G /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.genes.gff3 -o ../stringtie_merged.gtf list_to_merge.txt
                     conda deactivate"

 elif [ $1 -eq 2 ]; then
      sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$H \
      -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
   	--wrap ". $CONDA; conda activate newrnapipeline
                     gffcompare -r /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.genes.gff3 -o ../compared stringtie_merged.gtf
                     conda deactivate"

elif [ $1 -eq 3 ]; then
      sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$H \
      -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
   	--wrap ". $CONDA; conda activate newrnapipeline
                     chmod +x prepDE.py
                     module load python/2.7
                     ./prepDE.py -g ../gene_count_matrix.csv -i ./sample_list.txt
                     conda deactivate"
                     
                     
fi