
#!/bin/sh

F="/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads"
G="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"
##################################################################################################

if [ $1 -eq 1 ]; then
     sbatch --mem=300000 -N1 -n20 --ntasks-per-node=20 --workdir=$F --job-name "HaplotypeCaller" -o "$F/HaplotypeCaller.out" -e "$F/HaplotypeCaller.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "cd $F
             gatk HaplotypeCaller --reference $G --input R16b.SplitNCigarReads.split.bam --output R16b.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF
cd $F
             gatk HaplotypeCaller --reference $G --input T_M4C.SplitNCigarReads.split.bam --output T_M4C.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_M7.SplitNCigarReads.split.bam --output T_M7.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_R9D.SplitNCigarReads.split.bam --output T_R9D.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_M2.SplitNCigarReads.split.bam --output T_M2.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_M5.SplitNCigarReads.split.bam --output T_M5.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_R11B.SplitNCigarReads.split.bam --output T_R11B.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd $F
             gatk HaplotypeCaller --reference $G --input T_R15B.SplitNCigarReads.split.bam --output T_R15B.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 

             "

fi    