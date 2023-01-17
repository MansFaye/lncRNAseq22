#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --mail-user=mohamed.faye@students.unibe.ch
#SBATCH --job-name="Assembly_MF"
#SBATCH -c 2
#SBATCH --time=0:06:00
#SBATCH --mem=4G

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

if [ 1 == 0 ]
then
	for bam in $mapping_dir/*sorted*.bam
	do
		fname=$(basename $bam)
		ID=${fname:0:3}
		stringtie $bam -o $assembly_dir/${ID}Annot.gtf -p 6 -G $ref_annotation --rf 
	done
fi

# Merge the GTF files

if [ 1 == 0 ]
then
	stringtie --rf --merge $assembly_dir/*.gtf -o $assembly_dir/merged.gtf -G $ref_annotation
fi


# Map gene names to transcript names

if [ 1 == 0 ]
then
        awk -F $'\t' '$3=="transcript"' $assembly_dir/merged.gtf | grep 'gene_name' | awk -F $'\t' '{print $9}' | awk -F ';' '{print $2, $3}' | awk 'BEGIN { OFS="\t" }{print $2, $4}' | sed 's/["]//g' > $assembly_dir/IDs_map.txt
fi


# Summarize the content of the merged GTF gile

if [ 1 == 0 ]
then
	## Number of transcripts, genes, exons and novel transcripts
	echo "transcripts" | tr '\n' '\t' > $assembly_dir/summary.txt
	awk -F $'\t' '$3=="transcript"' $assembly_dir/merged.gtf | wc -l >> $assembly_dir/summary.txt

	echo "genes" | tr '\n' '\t' >> $assembly_dir/summary.txt
	awk -F $'\t' '$3=="transcript"{print $9}' $assembly_dir/merged.gtf | awk -F ';' '{print $1}' | sort | uniq -c | wc -l >> $assembly_dir/summary.txt

	echo "novel_transcripts" | tr '\n' '\t' >> $assembly_dir/summary.txt
	awk -F $'\t' '$3=="transcript"' $assembly_dir/merged.gtf | grep -v gene_name | wc -l >> $assembly_dir/summary.txt

	echo "exons" | tr '\n' '\t' >> $assembly_dir/summary.txt
	awk -F $'\t' '$3=="exon"' $assembly_dir/merged.gtf | wc -l >> $assembly_dir/summary.txt

	## Extract a list of single-exon transcripts 
	grep 'exon_number "2"' $assembly_dir/merged.gtf | awk -F $'\t' '{print $9}' | awk -F ';' '{print $2}' > $assembly_dir/multi_exon_tx.txt
	awk -F $'\t' '$3=="transcript"{print $9}' $assembly_dir/merged.gtf | awk -F ';' '{print $2}' > $assembly_dir/all_tx.txt
	diff --new-line-format="" --unchanged-line-format="" <(sort $assembly_dir/all_tx.txt) <(sort $assembly_dir/multi_exon_tx.txt) > $assembly_dir/single_exon_tx.txt
	rm $assembly_dir/multi_exon_tx.txt $assembly_dir/all_tx.txt

	## Add the # of single-exon transcripts to the summary file
	echo "single_exon_transcripts" | tr '\n' '\t' >> $assembly_dir/summary.txt
	cat $assembly_dir/single_exon_tx.txt | wc -l >> $assembly_dir/summary.txt

	echo "single_exon_novel_transcripts" | tr '\n' '\t' >> $assembly_dir/summary.txt
	grep MSTRG $assembly_dir/single_exon_tx.txt | wc -l >> $assembly_dir/summary.txt
fi


