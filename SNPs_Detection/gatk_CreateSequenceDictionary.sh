
O="/data/home/databases/Plutea_genome_reefgenomics"
R="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"

########################

if [ $1 -eq 1 ]; then
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$O --job-name "CreateSeqDict" -o "$O/CreateSequenceDictionary.out" -e "$O/CreateSequenceDictionary.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "gatk CreateSequenceDictionary -REFERENCE $R -OUTPUT ${R%*.*}.dict" 

elif [ $1 -eq 2 ]; then
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$O --job-name "CreateSeqDict" -o "$O/CreateSequenceDictionary.out" -e "$O/CreateSequenceDictionary.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "samtools faidx $R" 
						
fi                        