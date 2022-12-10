H="/data/home/Plutea_mangroves/output_Plutea_host/filtered/fastqc"

########################

if [ $1 -eq 1 ]; then
   mkdir -p $H
   sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$H \
   -p hive1d,hive7d,preempt7d,preempt31d,queen \
	--wrap ". $CONDA; conda activate multiqc_env; cd \"$H\"; multiqc .; conda deactivate"

fi