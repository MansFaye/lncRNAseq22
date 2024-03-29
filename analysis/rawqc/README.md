Quality Control of the raw reads is ran here using `doQC1.sh`. 
### Modules used
* Fastqc ver. 0.11.9
* MultiQC ver. 1.8
* Fastp (facultative)

### Input data
* RNA-seq data from TruSeq Stranded mRNA libraries, 3 replicates per clone.

### doQC1 script
The script first creates links to the sequencing data, making it more easily accessible.
It then runs the fastqc, storing results in the `./results` directory. Multiqc allows us to merge all results and view them in a single html file, which facilitates comparisons.
If the reads are of bad quality (for more info on fastqc's quality measures, see: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/), fastp can be used to filter reads or to trim eventual adapter sequences.

### Outputs, next steps
If the quality of reads is deemed acceptable, one can move on to mapping the reads to a reference genome.
MultiQC's per-base quality graph and read counts table are good figures to export for this step.
