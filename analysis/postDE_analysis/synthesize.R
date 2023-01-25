# Load Modules

library(ggplot2)
library(plyr)

setwd("/data/courses/rnaseq_course/lncRNAs/Project1/users/mfaye")

# Load all data 

novels_table <- read.csv("./analysis/Diff_Expression/DE_novel_transcripts.csv")

polya_intersects <- read.csv("./analysis/cpat/polya_transcript_ids.txt", sep =" ", header=F)
cage_intersects <- read.csv("./analysis/cpat/cage_transcript_ids.txt", sep =" ", header=F)
single_exon_tx <- read.csv("./analysis/cpat/single_exon_tx.txt", sep =" ", header=F)
pcp <- read.csv("./analysis/cpat/output", sep ="\t")

# Sifting through pcp results

pcp_noncoding <- pcp[which(pcp$coding_prob < 0.364),]

total_tx <- length(pcp$mRNA_size) # for the figure later
novel_tx <- length(grep("MSTRG", row.names(pcp), value = TRUE))
noncoding_tx <- length(pcp_noncoding$mRNA_size)
noncoding_novel_tx <- length(grep("MSTRG", row.names(pcp_noncoding), value = TRUE))

rm(pcp)

# Summarize all data in one place

results_tbl <- novels_table[,c(1,4,5)]
results_tbl$polya <- ifelse(results_tbl$target_id %in% polya_intersects$V3, TRUE ,FALSE)
results_tbl$cage <- ifelse(results_tbl$target_id %in% cage_intersects$V3, TRUE ,FALSE)
results_tbl$single_exon <- ifelse(results_tbl$target_id %in% single_exon_tx$V3, TRUE ,FALSE)
results_tbl$coding <- ifelse(results_tbl$target_id %in% row.names(pcp_noncoding), FALSE ,TRUE)

candidates <- results_tbl[which(!results_tbl$coding),]
top_candidates <- results_tbl[which(!results_tbl$coding & (results_tbl$cage | results_tbl$polya)),]

# Save the tables
write.csv(candidates, "./analysis/postDE_analysis/novel_candidates.csv", row.names=FALSE)
write.csv(top_candidates, "./analysis/postDE_analysis/top_novel_candidates.csv", row.names=FALSE)

# Make a df for plotting
df1 <- data.frame(coding_status=rep(c("coding", "non-coding"), each=2), annotation=c("novel","annotated","novel","annotated"),
                  num_tx=c(novel_tx - noncoding_novel_tx, (total_tx - novel_tx)-(noncoding_tx - noncoding_novel_tx), noncoding_novel_tx, noncoding_tx - noncoding_novel_tx))
df_sorted <- arrange(df1, annotation, num_tx) 
df_cumsum <- ddply(df_sorted, "annotation",
                   transform, label_ypos=cumsum(num_tx))

# Plot some general info

jpeg("./analysis/postDE_analysis/Protein_Coding.png", quality = 100)
ggplot(data=df_cumsum, aes(x=annotation, y=num_tx, fill=coding_status)) + 
  labs(title="Proportion of protein-coding transcripts", x="", y = "Number of transcripts") +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=num_tx), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()
dev.off()

