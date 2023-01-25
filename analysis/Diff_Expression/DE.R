# Load Modules

install.packages("BiocManager")
BiocManager::install("rhdf5", force=TRUE)
install.packages("devtools") 
library("devtools")
devtools::install_github("pachterlab/sleuth")
library(sleuth)
library(rhdf5)

setwd("/data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye")

# Define directory paths

base_dir <- "./data/quantif"
output_dir <- "./analysis/Diff_Expression"


# Load gene name/transcript ID mapping

ids <- read.table('./data/assembly/IDs_map.txt', sep="\t", header=FALSE)
names(ids) = c("target_id", "gene_name")

# Describe metadata

sample_id <- dir(file.path(base_dir))
kallisto_dirs <- file.path(base_dir, sample_id)

sample <- c("holoclonal1", "holoclonal2", "holoclonal5", "parental1", "parental2", "parental3")
condition <- c("H", "H", "H", "P", "P", "P")
s2c = data.frame(sample,condition)
s2c <- dplyr::mutate(s2c, path = kallisto_dirs)
print(s2c)

# Run sleuth

sleuth_out <- sleuth_prep(s2c, ~condition, target_mapping = ids, aggregation_column = 'gene_name', extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5), num_cores = 1)
## transcript-level only, for the volcano plot which doesn't support pval aggregation. 
#sleuth_out <- sleuth_prep(s2c, ~condition, target_mapping = ids, extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5), num_cores = 1)

sleuth_out <- sleuth_fit(sleuth_out, ~condition)

sleuth_out <- sleuth_wt(sleuth_out, 'conditionP') #approximates log2fold change

models(sleuth_out)

# Summarize results for significant transcripts + genes

sleuth_table_gn <- sleuth_results(sleuth_out, 'conditionP', 'wt', show_all = FALSE) #gene-level mode, aggregated pvalues
sleuth_table_tx <- sleuth_results(sleuth_out, 'conditionP', 'wt', show_all = FALSE, pval_aggregate = FALSE) #transcript level, observed counts(not tpm)

sleuth_table_gn <- dplyr::filter(sleuth_table_gn, qval <= 0.05)
sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
novels_table <- sleuth_table_tx[is.na(sleuth_table_tx$gene_name),] 

# Save tables

write.csv(sleuth_table_gn, "./analysis/Diff_Expression/DE_genes.csv", row.names=FALSE)
write.csv(sleuth_table_tx, "./analysis/Diff_Expression/DE_transcripts.csv", row.names=FALSE)
write.csv(novels_table, "./analysis/Diff_Expression/DE_novel_transcripts.csv", row.names=FALSE)

# Visualization
## Volcano plot
graph_points <- sleuth_results(sleuth_out, 'conditionP', 'wt', show_all = FALSE) # Make sure sleuth_prep was run without pval aggregation

volc_plot <- plot_volcano(sleuth_out, "conditionP", test_type = "wt", sig_level = 0.05, point_alpha = 0.2, sig_color = "red") +
  theme_bw() +
  geom_text(aes(graph_points$b, -log10(graph_points$qval), label=ifelse(qval < 1e-75, graph_points$gene_name,''), 
            hjust = 0, 
            vjust = 0,
            size=2), show.legend=FALSE)
dev.print(file="./analysis/Diff_Expression/volcano.png", device=png, width=800)

## Top candidates expression levels (see after step 7)
if (1==0){
  jpeg("./analysis/Diff_Expression/MSTRG.12388.2.png", quality = 100)
  plot_bootstrap(sleuth_out, 
                 target_id = "MSTRG.12388.2", 
                 units = "tpm", 
                 color_by = "condition")
  dev.off()
  
  jpeg("./analysis/Diff_Expression/MSTRG.10368.2.png", quality = 100)
  plot_bootstrap(sleuth_out, target_id = "MSTRG.10368.2", units = "tpm", color_by = "condition")
  dev.off()
  
  jpeg("./analysis/Diff_Expression/MSTRG.32577.2.png", quality = 100)
  plot_bootstrap(sleuth_out, target_id = "MSTRG.32577.2", units = "tpm", color_by = "condition")
  dev.off()
}


