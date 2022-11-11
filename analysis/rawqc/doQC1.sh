#!/bin/bash
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="FastQC"
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --mem=25G

# Move to project directory and set directory variables
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye

raw_reads_dir="../../../fastq"
holoclonal_reads=$raw_reads_dir/1*
parental_reads=$raw_reads_dir/P*
links_dir="./analysis/rawqc/fastqlnk"

# Load modules

module add UHTS/Quality_control/fastqc/0.11.9
module add UHTS/Analysis/MultiQC/1.8

# Create links to our fastqs in the wd and get read count for each replicate

if [ 1 == 0 ]
then
	for fastq in $holoclonal_reads
	do
        	echo $fastq >> ./analysis/rawqc/read_count.txt
        	ln -s ../../../$fastq $links_dir
        	zcat $fastq | awk '/@/' | wc -l >> ./analysis/rawqc/read_count.txt
	done

	for fastq in $parental_reads
	do
	        echo $fastq >> ./analysis/rawqc/read_count.txt
        	ln -s ../../../$fastq $links_dir
        	zcat $fastq | awk '/@/' | wc -l >> ./analysis/rawqc/read_count.txt
	done
fi 

# Run the fastQC

if [ 1 == 0 ]
then
	mkdir ./analysis/rawqc/results
	fastqc -t 4 -o ./analysis/rawqc/results ./analysis/rawqc/fastqlnk/* 
	#fastp -i bc2023.fastq.gz -W 10 -5 -3 -M 10 -l 1000 -G -A -w 4 -o bc2023_clean.fastq.gz
fi

# Run multiQC

if [ 1 == 0 ]
then
        cd ./analysis/rawqc/results
        multiqc .
fi








