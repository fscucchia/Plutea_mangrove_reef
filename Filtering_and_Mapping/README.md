
The following document contains the bioinformatic pipeline used for cleaning, aligning and assembling _P. lutea_ raw RNA sequences.

---

**Tools used**  

Quality check: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)  
Quality trimming: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
Alignment to the reference genome: [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)  
Preparation of alignment for the assembly: [SAMtools](http://www.htslib.org/doc/samtools.html)  
Transcripts assembly and quantification: [StringTie](https://ccb.jhu.edu/software/stringtie/) 

---

#### Create and activate a new conda environment

```
conda create -n newrnapipeline
conda activate newrnapipeline
```

#### Install all necessary programs within your new conda environment

```
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda multiqc  
conda install -c bioconda/label/broken trimmomatic
conda install -c bioconda cutadapt
conda install -c bioconda fastqc 

wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.5.Linux_x86_64.tar.gz
tar xzf stringtie-2.1.5.Linux_x86_64.tar.gz
```

### Quality filtering 

[Design table](https://github.com/fscucchia/Plutea_mangrove_reef/blob/main/Metadata/design_PE.csv) with samples list.

Run [`fastq-filter-PE_Conda_1_RUN1.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/fastq-filter-PE_Conda_1_RUN1.sh) that performs quality filtering. Argument 1 uses FastQC to perform the initial quality check of raw reads. Argument 2 calls for the script [`fastq-filter_Conda_job_1.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/fastq-filter_Conda_job_1.sh), which contains all the cutadapt and trimmomatics commands, and removes the adapters using the file [`adapters.fa`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/adapters.fa). Argument 3 uses FastQC to perform the quality check of the filtered reads.
Run [`MultiQC.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/MultiQC.sh) for multiQC.

### Alignment of the clean reads to _P. lutea_ reference genome 

Create a new directory for Hisat.
```
mkdir HISAT
```
Run the script [`StyHISAT2.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/StyHISAT2.sh) that does the following:
1) Symbolically links the filtered fastq files to the HISAT directory.
```
/data/home/Plutea_mangroves/output_Plutea_host/fastqc_filtered/*gz.filtered ./
```
2) Index the _P. lutea_ reference genome in the reference directory and performs the alignment using the script [`StyHISAT_withSummary.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/StyHISAT_withSummary.sh). 
```
hisat2-build -f /data/home/databases/Plutea_genome_reefgenomics/plut_final_2.1.fasta ./Plut_ref
```
```
#StyHISAT_withSummary.sh:

   F="/data/home/Plutea_mangroves/output_Plutea_host/filtered"
   # Aligning paired end reads
   # The R1 in array1 is changed to R3 in the for loop. SAM files are of both forward and reverse reads
   array1=($(ls $F/*_R1.fastq.gz.filtered))

   # Bam files are created and sorted, since Stringtie takes sorted file as input
   # The sam file is removed at the end since it is not needed anymore
   # The command --summary-file ${i}.txt reates a summary file per sample, which can be used by multiqc

   for i in ${array1[@]}; do
        hisat2 -p 8 --new-summary --rf --dta -q -x Plut_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R3/) -S ${i}.sam --summary-file ${i}.txt 
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
   done
```
3) Runs multiqc at the end of mapping if the `--summary-file` option was enabled. 

<details>
<summary>Troubleshooting tips!</summary>
<br>
After the alignment of the first sample, I got this message "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory"
So I ran these 3 commands: $ conda config --add channels bioconda
                           $ conda config --add channels conda-forge
                           $ conda install samtools==1.11
</details>

### Assemble aligned reads and quantify transcripts 

Create a new directory for StringTie. 
```
mkdir StringTie
```
Run the script [`StySTRINGTIE.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/StySTRINGTIE.sh)
that does the following:
1) Creates a symbolic link to the reference genome gff file inside the stringtie directory
```
ln -s /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.genes.gff3 ./
```
2) Runs the script [`StyStringTie_assembly.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/blob/main/Filtering_and_Mapping/StyStringTie_assembly.sh) to assemble the aligned reads and quantify transcripts:
```
  ##!/bin/bash

  #Specify working directory
  W="/data/home/Plutea_mangroves/output_Plutea_host/filtered"

  #StringTie reference-guided assembly
  #These BAM files contain both forward and reverse reads
  array1=($(ls $W/*_R1.gz.filtered.bam))

  for i in ${array1[@]}; do
        stringtie -A gene_abundance/{i}.gene_abund.tab -p 8 --rf -e -G Pastreoides_all_v1.gff -o ${i}.gtf ${i}
        mv ${i}.gtf /data/home/Plutea_mangroves/output_Plutea_host/StringTie/BAM
        echo "StringTie-assembly-to-ref ${i}" $(date)
  done
```

### Assess the performance of the assembly 

Install gffcompare within your conda environment.
```
conda activate newrnapipeline
conda install -c bioconda gffcompare
```
Then create a .txt file ([`list_to_merge.txt`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/list_to_merge.txt)) listing all of the file names to be merged. This file needs to be in the StringTie directory.

Run the script [`Stringtie_merge_compare_count.sh`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/Stringtie_merge_compare_count.sh), which does the following:
1) Merges the GTF files generated from the assembly to assess how well the predicted transcripts track to the reference annotation gff file.
```
stringtie --merge -p 8 -G /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.genes.gff3 -o ../stringtie_merged.gtf list_to_merge.txt
```
2) Uses the program gffcompare to compare the merged GTF files to the reference genome.
```
gffcompare -r /data/home/databases/Plutea_genome_reefgenomics/plut2v1.1.genes.gff3 -o ../compared stringtie_merged.gtf
```
3) Compiles the GTF files into gene and transcript count matrices. The StringTie program includes the script `prepDE.py` that compiles the assembly-generated files into gene and transcript count matrices. This script requires as input a list of sample names and their file paths which has to be manually created. This .txt file ([`sample_list.txt.`](https://github.com/fscucchia/Plutea_mangrove_reef/tree/main/Filtering_and_Mapping/sample_list.txt) needs to be in the StringIie directory.
```
./prepDE.py -g ../gene_count_matrix.csv -i ./sample_list.txt
```
<details>
<summary>Troubleshooting tips!</summary>
<br>
I got syntax related-errors when running the ./prepDE.py with python3. So I used '$ module load python/2.7' to run the script.                     
</details>
