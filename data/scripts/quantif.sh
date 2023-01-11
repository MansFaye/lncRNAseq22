#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="Quantif_MF"
#SBATCH -c 8
#SBATCH --time=4:30:00
#SBATCH --mem=60G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye

assembly_dir="./data/assembly"
quantification_dir="./data/quantif"
fastq_links_dir="./analysis/rawqc/fastqlnk"
ref_files_dir="../../references"

mkdir $quantification_dir

# Load modules
module add UHTS/Analysis/kallisto/0.46.0
module load UHTS/Assembler/cufflinks/2.2.1

# Generate a FASTA from our GTF and ref. genome

if [ 1 == 1 ]
then
	gffread -w $quantification_dir/transcripts.fa -g $ref_files_dir/GRCh38.genome.fa $assembly_dir/merged.gtf
fi

# Indexing the FASTA file

if [ 1 == 1 ]
then
	kallisto index -i $quantification_dir/transcripts_kallisto $quantification_dir/transcripts.fa
fi 

# Quantification
## Parental reads

if [ 1 == 1 ]
then
	for i in {1..3}
	do
		kallisto quant -i $quantification_dir/transcripts_kallisto -b 100 -o $quantification_dir/parental${i} -t 8 --rf-stranded $fastq_links_dir/P${i}*R1*.fastq.gz $fastq_links_dir/P${i}*R2*.fastq.gz 
	done
fi

## Holoclonal reads

if [ 1 == 1 ]
then
        for i in {1,2,5}
        do
                kallisto quant -i $quantification_dir/transcripts_kallisto -b 100 -o $quantification_dir/holoclonal${i} -t 8 --rf-stranded $fastq_links_dir/1_${i}*R1*.fastq.gz $fastq_links_dir/1_${i}*R2*.fastq.gz
        done
fi


