
S="/data/home/Plutea_mangroves/scripts"
F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"

########################

if [ $1 -eq 1 ]; then
     mkdir -p $S
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$S --job-name "MergeBamAlignment" -o "$F/MergeBamAlignment3.out" -e "$F/MergeBamAlignment3.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "chmod u+x MergeBamAlignment.sh
			 ./MergeBamAlignment.sh" 
					  						
fi                        