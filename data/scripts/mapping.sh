#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="STAR_mapping"
#SBATCH -c 6
#SBATCH --time=2:00:00
#SBATCH --mem=150G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye # Move to project directory

rawqc_links_dir="./analysis/rawqc/fastqlnk"
ref_files_dir="../../references"
mapping_dir="./data/mapping"
indexes_dir="$ref_files_dir/indexes"
forward_fastqs=$rawqc_links_dir/*R2*
mkdir $indexes_dir

#Load Modules
module add UHTS/Aligner/STAR/2.7.9a
module add UHTS/Analysis/samtools/1.8

# Index the reference genome using STAR

if [ 1 == 0 ]
then
	STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $indexes_dir --genomeFastaFiles $ref_files_dir/GRCh38.genome.fa 
fi

# Align the reads using STAR

if [ 1 == 1 ]
then
	for fastq in $forward_fastqs
	do
		fname=$(basename $fastq .gz)
		ID=${fname:0:3}
		reverse=$rawqc_links_dir/$ID*R1*
		STAR --runMode alignReads --genomeDir $indexes_dir --runThreadN 6 --genomeLoad  LoadAndKeep --limitBAMsortRAM 30000000000 --readFilesCommand zcat --readFilesIn $fastq $reverse --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $mapping_dir/$ID
	done
fi

# Index one of the BAM files to visualize it in IGV 

if [ 1 == 1 ]
then
	bam=$mapping_dir/1_1*sorted*.bam # BAM file to be indexed
	samtools index -b -@ 6 $bam
fi





