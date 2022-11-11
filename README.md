# Annotation and characterization of lncRNAs in Lung Cancer
Identification of differentially expressed protein-coding and long noncoding genes in the cancer line A549's parental and
holoclonal lines.  
-insert some rendition of the abstract-
See also: 1.Tièche, C. C. et al. Tumor Initiation Capacity and Therapy Resistance Are Differential Features of EMT-Related Subpopuliations in the NSCLC Cell Line A549. Neoplasia 21, 185-196 (2019). doi:10.1016/j.neo.2018.09.008

### Main Steps:
1. Quality Control of raw reads
2. Mapping
3. Transcriptome assembly
4. Quantification
5. Differential Expression
6. Integrative analysis & Prioritization

## Packages
FastQC (ver. 0.11.9):
MultiQC (ver. 1.8):

## Data
### Experimental data:
RNA-seq data from TruSeq Stranded mRNA libraries, 3 replicates per clone.

### Reference data:


# Quality Control of raw reads
We check for quality of the sequencing data using the FastQC package, and visualize results for all replicates using MultiQC.
See README in `./analysis/rawqc` for detailed information.

# Read Mapping to the ... genome

