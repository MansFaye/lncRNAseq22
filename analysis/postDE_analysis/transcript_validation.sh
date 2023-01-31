#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="Beds_MF"
#SBATCH -c 3
#SBATCH --time=1:00:00
#SBATCH --mem=30G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye # Move to project directory

beds_dir="./analysis/postDE_analysis/bed_files"
assembly_dir="./data/assembly"
annotation=$assembly_dir/merged.gtf
output_dir="./analysis/postDE_analysis"
cpat_dir=$output_dir/cpat
quantif_dir="./data/quantif"

cage_bed=$beds_dir/hg38_fair+new_CAGE_peaks_phase1and2.bed
polya_bed=$beds_dir/atlas.clusters.2.0.GRCh38.96.bed

mkdir $cpat_dir

# Load Modules

module load UHTS/Analysis/BEDTools/2.29.2
module load SequenceAnalysis/GenePrediction/cpat/1.2.4
module load R/latest

# Make a BED file of transcript start sites in the GTF annotation file from step 3.

if [ 1 == 1 ]
then
	awk -F $'\t' '$3=="transcript"' $assembly_dir/merged.gtf | awk -F $'\t' 'BEGIN { OFS="\t" }{print $1, $4, $5, $9, ".", $7}' | awk -F $'\t' 'BEGIN { OFS="\t" }{if ($6 == "-") {print $1, $3-50, $3+50, $4, $5, $6} else {print $1, $2-50, $2+50, $4, $5, $6}}' > $beds_dir/annotation_start.bed
fi

# Make a BED file of transcript end sites in the GTF annotation file from step 3.

if [ 1 == 1 ]
then
        awk -F $'\t' '$3=="transcript"' $assembly_dir/merged.gtf | awk -F $'\t' 'BEGIN { OFS="\t" }{print $1, $4, $5, $9, ".", $7}' | awk -F $'\t' 'BEGIN { OFS="\t" }{if ($6 == "-") {print $1, $2-50, $2+50, $4, $5, $6, ".", ".", ".", ".", ".", "."} else {print $1, $3-50, $3+50, $4, $5, $6, ".", ".", ".", ".", ".", "."}}' | grep '^chr' > $beds_dir/annotation_end.bed
fi

# Find overlap between CAGE clusters and our start site BED file.
if [ 1 == 1 ]
then
	bedtools intersect -s -wa -a $beds_dir/annotation_start.bed -b $cage_bed > $beds_dir/cage_intersects.bed
fi

# Find overlap between PolyA Site clusters and our end site BED file.
if [ 1 == 1 ]
then
	sed -e 's/^/chr/' $polya_bed > $beds_dir/polya_clusters_corrected.bed # Chr names from this BED file dont start with chr.
	bedtools intersect -s -wa -a $beds_dir/annotation_end.bed -b $beds_dir/polya_clusters_corrected.bed > $beds_dir/polyA_intersects.bed
fi

# Extract transcript IDs of all novel transcripts that intersect with CAGE peaks and PolyA clusters
if [ 1 == 1 ]
then
	awk -F $'\t' '{print $4}' $beds_dir/cage_intersects.bed | grep -v gene_name | awk -F ';' '{print $2}' | sort | uniq > $output_dir/cage_transcript_ids.txt
	awk -F $'\t' '{print $4}' $beds_dir/polyA_intersects.bed | grep -v gene_name | awk -F ';' '{print $2}' | sort | uniq > $output_dir/polya_transcript_ids.txt
fi

# Evaluate Protein-coding potential (ORFs)
if [ 1 == 1 ]
then
	cpat.py -x $cpat_dir/Human_Hexamer.tsv -d $cpat_dir/Human_logitModel.RData -g $quantif_dir/transcripts.fa -o $cpat_dir/output
fi


