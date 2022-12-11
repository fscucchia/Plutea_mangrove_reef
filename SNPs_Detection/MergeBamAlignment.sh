
#!/bin/sh

# F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
# G="/data/home/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"
#   array1=($(ls $F/*_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam))
#   for i in ${array1[@]}; do
#                 gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM ${i}.FastqToSam.unmapped.rg.bam --ALIGNED_BAM ${i}_R1_concat.fastq.gz.filtered.bam --OUTPUT ${i}.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
# 				touch ${i}.MergeBamAlignment.done

#   done

G="/data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta"


cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered_single.bam --OUTPUT R16b.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch R16b.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_M2.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_M2.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_M4C.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_M4C.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_M5.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_M5.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_M7.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_M7.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_R11B.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_R11B.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_R15B.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_R15B.MergeBamAlignment.done

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM /data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered_single.bam --OUTPUT T_R9D.fastq.gz.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch T_R9D.MergeBamAlignment.done
