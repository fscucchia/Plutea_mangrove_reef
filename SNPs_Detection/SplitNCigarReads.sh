
#!/bin/sh

G="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
O="/data/home/Plutea_mangroves/output_Plutea_host/SplitNCigarReads/all_rRNA_cleaned"
R="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=300000 --workdir=$G --job-name "SplitNCigarReads" -o "$G/SplitNCigarReads.out" -e "$G/SplitNCigarReads.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
     --wrap "cd $G
             gatk SplitNCigarReads -R $R -I R16b_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/R16b.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_M4C_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_M4C.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_M7_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_M7.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_R9D_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_R9D.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_M2_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_M2.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_M5_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_M5.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_R11B_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_R11B.SplitNCigarReads.split.bam
cd $G
             gatk SplitNCigarReads -R $R -I T_R15B_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/T_R15B.SplitNCigarReads.split.bam"

fi 
