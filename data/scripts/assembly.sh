#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="Assembly_MF"
#SBATCH -c 6
#SBATCH --time=2:30:00
#SBATCH --mem=100G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye # Move to project directory

ref_files_dir="../../references"
mapping_dir="./data/mapping"
assembly_dir="./data/assembly"
ref_annotation=$ref_files_dir/gencode.v21.chr_patch_hapl_scaff.annotation.gtf

mkdir $assembly_dir

# Load Modules

module add UHTS/Aligner/stringtie/1.3.3b
module add UHTS/Analysis/samtools/1.8

# Perform the assembly

if [ 1 == 1 ]
then
	for bam in $mapping_dir/*sorted*.bam
	do
		fname=$(basename $bam)
		ID=${fname:0:3}
		stringtie $bam -o $assembly_dir/${ID}Annot.gtf -p 6 -G $ref_annotation --rf 
	done
fi

# Merge the GTF files

if [ 1 == 1 ]
then
	stringtie --rf --merge $assembly_dir/*.gtf -o $assembly_dir/merged.gtf -G $ref_annotation
fi


