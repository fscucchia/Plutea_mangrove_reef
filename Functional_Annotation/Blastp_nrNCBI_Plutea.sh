
H="/data/home/databases/nr_NCBI"
Q="/data/home/databases/Plutea_genome_reefgenomics"
O="/data/home/Plutea_mangroves/output_Plutea_host/functional_annot_Blastp"
#########################

if [ $1 -eq 1 ]; then
   sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --job-name blastpPlutea --workdir=$H \
   -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
 	--wrap ". $CONDA; conda activate functional_annotation; gunzip -c nr.gz; conda deactivate"

elif [ $1 -eq 2 ]; then
   sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --job-name blastpPlutea --workdir=$H \
   -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
 	--wrap ". $CONDA; conda activate functional_annotation; makeblastdb -in nr.fasta -dbtype prot; conda deactivate"

elif [ $1 -eq 3 ]; then
   sbatch --cpus-per-task=30 --mem=500000 --job-name blastpPlutea -o "$O/blastp2.out" -e "$O/blastp2.err" --workdir=$H \
   -p queen \
--wrap ". $CONDA; conda activate functional_annotation; blastp -query $Q/plut2v1.1.proteins.fasta -db $H/nr.fasta -num_threads 20 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out $O/Plut.annot_blastp.23012022; conda deactivate"

fi


