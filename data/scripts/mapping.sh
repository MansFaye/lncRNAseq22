#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="STAR_mapping"
#SBATCH -c 6
#SBATCH --time=5:00:00
#SBATCH --mem=100G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye

rawqc_links_dir="./analysis/rawqc/fastqlnk"
ref_files_dir="../../references"
mapping_dir="./data/mapping"

mkdir $ref_files_dir/indexes

indexes_dir="$ref_files_dir/indexes"

#Load Modules
module add UHTS/Aligner/STAR/2.7.9a

# Index the reference genome using STAR

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $indexes_dir \
--genomeFastaFiles $ref_files_dir/GRCh38.genome.fa \

