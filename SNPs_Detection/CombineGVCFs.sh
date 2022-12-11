
#!/bin/sh

N="/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads"
REF="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=300000 --workdir=$N --job-name "CombineGVCFs" -o "$N/CombineGVCFs.out" -e "$N/CombineGVCFs.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk CombineGVCFs -R "${REF}" -V R16b.HaplotypeCaller.g.vcf.gz -V T_M2.HaplotypeCaller.g.vcf.gz -V T_M5.HaplotypeCaller.g.vcf.gz -V T_M7.HaplotypeCaller.g.vcf.gz -V T_M4C.HaplotypeCaller.g.vcf.gz -V T_R11B.HaplotypeCaller.g.vcf.gz -V T_R15B.HaplotypeCaller.g.vcf.gz -V T_R9D.HaplotypeCaller.g.vcf.gz -O $N/cohort.g.vcf.gz"
fi




 
