

# files:
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
title[0]="T_R11B"
for1[0]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R1.fastq.gz.filtered"
rev1[0]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R11B_HTKTTDSX3_AGCGCTGTGT-GTTCCGAACG_L004_R2.fastq.gz.filtered"
title[1]="T_M5"
for1[1]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R1.fastq.gz.filtered"
rev1[1]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M5_HTKTTDSX3_ACTTACTTCA-TAATTCCAGA_L004_R2.fastq.gz.filtered"
title[2]="T_R9D"
for1[2]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R1.fastq.gz.filtered"
rev1[2]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R9D_HTKTTDSX3_GATATCACAC-ACATATAGTC_L004_R2.fastq.gz.filtered"
title[3]="R16b"
for1[3]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R1.fastq.gz.filtered"
rev1[3]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/R16b_HTKTTDSX3_GATAGCCTTG-CATTGTGGTA_L004_R2.fastq.gz.filtered"
title[4]="T_M7"
for1[4]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R1.fastq.gz.filtered"
rev1[4]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M7_HTKTTDSX3_CGCATTCCGT-CTCCACTAAT_L004_R2.fastq.gz.filtered"
title[5]="T_M4C"
for1[5]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R1.fastq.gz.filtered"
rev1[5]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M4C_HTKTTDSX3_GCCTAACGTG-CTTGAGAGCT_L004_R2.fastq.gz.filtered"
title[6]="T_M2"
for1[6]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R1.fastq.gz.filtered"
rev1[6]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_M2_HTKTTDSX3_CTTAAGTCGA-ACAACTACTG_L004_R2.fastq.gz.filtered"
title[7]="T_R15B"
for1[7]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R1.fastq.gz.filtered"
rev1[7]="/data/home/Plutea_mangroves/output_Plutea_host/filtered/T_R15B_HTKTTDSX3_TCACCGCGCT-CTAGTATCGA_L004_R2.fastq.gz.filtered"


#run:
#for i in ${!title[@]}; do
#	echo $i
#	echo ${title[i]}
#	echo ${for1[i]}
#	echo ${rev1[i]}
#	echo '----------------'
#done
