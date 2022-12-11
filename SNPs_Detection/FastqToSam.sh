
#!/bin/sh

# F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
# aarray1=($(ls $F/*_R1.fastq.gz.filtered))
# for i in ${array1[@]}; do
#         gatk FastqToSam \
# 				--FASTQ ${i} \
# 				--FASTQ2 $(echo ${i}|sed s/_R1/_R3/) \
# 				--OUTPUT ${i}.FastqToSam.unmapped.bam \
# 				--SAMPLE_NAME ${i}; touch ${i}.FastqToSam.done
# done

F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
array1=($(ls $F/*_R1.fastq.gz.filtered))
for i in ${array1[@]}; do
        gatk FastqToSam \
				--FASTQ ${i} \
				--OUTPUT ${i}.FastqToSam.unmapped.bam \
				--SAMPLE_NAME ${i}; touch ${i}.FastqToSam.done

done