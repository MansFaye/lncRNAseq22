#!/bin/bash
#SBATCH --job-name="sift"
#SBATCH --cpus-per-task=2
#SBATCH --time=0:02:00
#SBATCH --mem=5G

# Move to project directory and set directory variables

cd /data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye

mapping_dir="./data/mapping"
mapping_stats_dir="./analysis/mapping_stats"

output_csv=$mapping_stats_dir/stats.csv
temporary=$mapping_stats_dir/temp.csv

# Extract the mapped reads counts from the Log files

for logfile in $mapping_dir/*Log.final*
do
	fname=$(basename $logfile)
	ID=${fname:0:3}
	awk -v logfile="$logfile" -F $'\t' 'NR == 6 || NR == 9 || NR == 24 {print $2}' $logfile | tr '\n' '\t' | awk -v clone="$ID" 'BEGIN{OFS="\t"}{print clone, $0}' >> $temporary
done

# We now have IDs, total reads, uniquely aligned reads and multiply aligned reads. We need to do a few modifications

awk -F $'\t' 'BEGIN{OFS="\t"}{print $1, $2, $3+$4, ($3+$4)/$2*100, $3, $3/$2*100}' $temporary | awk 'BEGIN{print "sample\ttotal_reads\ttotal_mapped_reads\tmapped_reads_percentage\tuniquely_mapped_reads\tuniquely_mapped_reads_percentage"}1' > $output_csv 
rm $temporary

 
