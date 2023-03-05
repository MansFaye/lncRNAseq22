# Annotation and characterization of lncRNAs in Lung Cancer
Lung cancer is the first cause of cancer-related mortality, and non-small cell lung cancers (NSCLC) make up more than 80% of them. Three distinct subpopulations were identified within the NSCLC cell line A549, with one in particular showing high expression of cancer stem-cell (CSC)-related mRNAs and protein markers. CSCs have been known to be resistant to most conventional therapeutics, making them promising targets to improve cancer therapy. As the implication of lncRNAs in cancer stemness and tumorigenesis is becoming clearer in recent years, we identified protein-coding and lncRNA genes that were differentially expressed between the holoclonal subpopulation and parental line. 

See also: Tièche, C. C. et al. Tumor Initiation Capacity and Therapy Resistance Are Differential Features of EMT-Related Subpopuliations in the NSCLC Cell Line A549. Neoplasia 21, 185-196 (2019). doi:10.1016/j.neo.2018.09.008

### Main Steps:
1. Quality Control of raw reads
2. Mapping
3. Transcriptome assembly
4. Quantification
5. Differential Expression
6. Integrative analysis & Prioritization

## Softwares / Packages
BEDTools (ver. 2.29.2): https://bedtools.readthedocs.io/en/latest/ 

CPAT (ver. 1.2.4): https://cpat.readthedocs.io/en/latest/

Cufflinks (ver. 2.2.1): http://cole-trapnell-lab.github.io/cufflinks/

FastQC (ver. 0.11.9): https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Integrative Genomics Viewer (ver. 2.15.2): https://software.broadinstitute.org/software/igv/

Kallisto (ver. 0.46.0): https://pachterlab.github.io/kallisto/

MultiQC (ver. 1.8): https://multiqc.info/

Samtools (ver. 1.8): https://www.htslib.org/

Sleuth (ver. 0.30.1): https://pachterlab.github.io/sleuth/

STAR (ver. 2.7.9a): https://github.com/alexdobin/STAR

Stringtie (ver. 1.3.3b): https://ccb.jhu.edu/software/stringtie/index.shtml

## Data
### Experimental data:
RNA-seq data from TruSeq Stranded mRNA libraries, 3 replicates per clone.

### Reference data:
Genome Reference Consortium Human Build 38, Fasta and annotation : Fetched from https://www.gencodegenes.org/human/release_21.html on 14.11.2022.
Note: The annotation used is the full comprehensive gene annotation, labeled as "ALL" on the source website.

FANTOM5 CAGE peaks, Human Build 38: Fetched from https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz on 03.01.2023.

PolyASite clusters v2.0, Human Build 38: Fetched from https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz on 03.01.2023.

CPAT Hexamer frequency table and logistic regression model: fetched from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/ on 12.01.2023.
 
# Quality Control of raw reads
We check for quality of the sequencing data using the FastQC package, and visualize results for all replicates using MultiQC.
See README in `analysis/rawqc` for detailed information.

# Mapping reads to the reference genome
We generate indexes for the reference genome and align the reads using STAR. 
Reminder: we have paired-end reads, which must be loaded together for alignement.
We then generate a bai index for at least one of our BAM files, in order to visualize the mapping in IGV.
See `data/scripts/mapping.sh` for details, and `data/mapping/` for the output files.

Before moving on, we check the proportion of mapped reads as well as uniquely aligned reads in STAR's Log files. There likely is a problem if the total number of aligned reads amounts to less than 70% of all the reads. This information is summarized using `analysis/mapping_stats/make_table.sh`


# Transcriptome Assembly
We create a reference-guided assembly using StringTie in order to include transcripts that are absent from the GENCODE annotation in our analysis. The assembly is performed once per replicate, and the results merged using StringTie's `--merge` function. We then extract a mapping of gene name to transcript ID for the Differential Expression step, as well as the number of genes, transcripts, exons, novel transcripts, and single-exon transcripts.
See `data/scripts/assembly.sh` for details, and `data/assembly/` for the output files


# Quantification
We first generate a FASTA file from our annotation file using the `gffread` function from Cufflinks. This will serve as a list of target sequences. Kallisto then indexes the FASTA file, before separately quantifying gene expression for each replicate. Kallisto outputs abundance data in plain text and HDF5 format, which can be found in `data/quantif/`.
See `data/scripts/quantif.sh` for details.


# Differential Expression
The Kallisto results are read by Sleuth, which compares transcript-level count estimates between holoclonal and parental replicates. It automatically applies the Benjamini-Hochberg procedure to control the false-discovery rate. Gene-level differential expression results are obtained by aggregating transcript-level results for genes that have multiple transcripts. Results for both genes and transcripts with q-value < 0.05 can be found in `analysis/Diff_Expression/`
See `analysis/Diff_Expression/DE.R` for details.


# Integrative Analysis
## Checking Transcription Start Site and Polyadenylation Site accuracy
We compare the transcription start sites from the previously generated assembly with peaks identified by CAP Analysis of Genome Expression (CAGE), using data from the FANTOM consortium (see Reference Data section). 
The same thing is done for the end of the transcripts (the Polyadenylation Sites), using the curated PolyASite clusters from the University of Basel (see Reference Data section).
The goal being to validate the novel (non-annotated) transcripts. We are more confident in the validity of transcripts whose start/end site overlap with CAGE peaks and/or PolyASite clusters.
See `analysis/postDE_analysis/transcript_validation.sh` for more details
## Checking protein-coding potential
In order to identify which of the novel transcripts are lncRNAs, we calculate protein-coding potential using cpat. This step is part of the `analysis/postDE_analysis/transcript_validation.sh` script.
The CPAT results, PolyA and CAGE intersects, and the differential expression results are all summarized using the `analysis/postDE_analysis/synthesize.R` R script.
The lncRNA candidates are listed in `analysis/postDE_analysis/novel_candidates.tsv`.


