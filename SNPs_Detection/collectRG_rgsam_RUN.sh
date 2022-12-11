
S="/data/home/Plutea_mangroves/scripts"
F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"

########################

if [ $1 -eq 1 ]; then
     mkdir -p $S
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$S --job-name "collectRG_rgsam" -o "$F/collectRG_rgsam.out" -e "$F/collectRG_rgsam.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "chmod u+x collectRG_rgsam.sh
			 ./collectRG_rgsam.sh" 
					  						
fi                        