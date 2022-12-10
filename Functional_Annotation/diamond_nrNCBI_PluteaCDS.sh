
O="/data/home/Plutea_mangroves/output_Plutea_host/functional_annot_DIAMOND"

########################à

if [ $1 -eq 1 ]; then
   sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --job-name diamondblast \
   -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
	--wrap ". $CONDA; conda activate functional_annotation; diamond blastp -d /data/home/databases/nr_NCBI/nr.dmnd -q /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.proteins.fasta -o $O/Plut.diamondBlastpNCBInr_20012022 -f 6 -b 20 --more-sensitive -e 0.00001 -k1; conda deactivate"

fi