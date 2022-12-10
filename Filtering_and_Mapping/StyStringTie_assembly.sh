#!/bin/bash

W="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
array1=($(ls $W/*_R1.gz.filtered.bam))
for i in ${array1[@]}; do
        stringtie -A gene_abundance/{i}.gene_abund.tab -p 8 --rf -e -G Pastreoides_all_v1.gff -o ${i}.gtf ${i}
        mv ${i}.gtf /data/home/Plutea_mangroves/output_Plutea_host/StringTie/BAM
        echo "StringTie-assembly-to-ref ${i}" $(date)
done