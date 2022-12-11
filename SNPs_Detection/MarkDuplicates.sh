
#!/bin/sh

G="/data/home/Plutea_mangroves/output_Plutea_host/filtered"

if [ $1 -eq 1 ]; then
     sbatch --mem=300000 -N1 -n20 --ntasks-per-node=20 --workdir=$G --job-name "MarkDuplicates" -o "$G/MarkDuplicates.out" -e "$G/MarkDuplicates.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT R16b.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT R16b_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE R16b_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_M4C.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_M4C_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_M4C_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_M7.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_M7_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_M7_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_R9D.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_R9D_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_R9D_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_M2.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_M2_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_M2_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_M5.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_M5_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_M5_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_R11B.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_R11B_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_R11B_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
             gatk MarkDuplicates --INPUT T_R15B.fastq.gz.MergeBamAlignment.merged.bam --OUTPUT T_R15B_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE T_R15B_MergeBamAlignment.MarkDuplicates.metrics"
fi   