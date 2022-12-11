
#!/bin/sh

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o R16b_rg.txt
samtools view -h R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r R16b_rg.txt |
  samtools view -b - > R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_M2_rg.txt
samtools view -h T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam| 
  rgsam tag -r T_M2_rg.txt |
  samtools view -b - > T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_R9D_rg.txt
samtools view -h T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_R9D_rg.txt |
  samtools view -b - > T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_M4C_rg.txt
samtools view -h T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_M4C_rg.txt |
  samtools view -b - > T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_M5_rg.txt
samtools view -h T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_M5_rg.txt |
  samtools view -b - > T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_M7_rg.txt
samtools view -h T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_M7_rg.txt |
  samtools view -b - > T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_R11B_rg.txt
samtools view -h T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_R11B_rg.txt |
  samtools view -b - > T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/Plutea_mangroves/output_Plutea_host/filtered
samtools view T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam -o T_R15B_rg.txt
samtools view -h T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r T_R15B_rg.txt |
  samtools view -b - > T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

