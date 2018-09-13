# scICW
scICW is an integrated pipeline for both gene expression and genomic architecture analysis in single cells. Our pipeline reveals the global expression profile of the populations, and also identifies the changes in transcriptome/genome including alternative splicing (AS), single-nucleotide polymorphisms (SNPs), RNA editing and gene fusion.  
# Dependencies
scICW requires the following softwares and reference data:  
● Samtools (https://github.com/samtools/samtools)  
● Vcftools (https://vcftools.github.io/index.html)  
● Fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
● Qualimap (http://qualimap.bioinfo.cipf.es/)  
● Multiqc (http://multiqc.info/)  
● Cutadapt (https://github.com/marcelm/cutadapt)  
● STAR (https://github.com/alexdobin/STAR)  
● RSEM (https://github.com/deweylab/RSEM)  
● Seurat (https://satijalab.org/seurat/)  
● Picard (https://broadinstitute.github.io/picard/)  
● GATK (https://software.broadinstitute.org/gatk/)  
● Sifit (https://bitbucket.org/hamimzafar/sifit)  
● StarFusion (https://github.com/STAR-Fusion/STAR-Fusion)  
● red_ML (https://github.com/BGIHENG/RED-ML)  

● QM_GTF ()  
● STAR_index ()  
● RSEM_index ()  
● Ref_fasta (hg38) ()  
● Ref_len ()  
● dbSNP ()   
● CTAT_lib ()  
● SimpleRepeat (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database)   
● ALU (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database)  
# Installation
scICW.v1.pl is provided in the ```bin``` directory.  
# Usage
## Parameters
```
Basic:
            * --conf     Configure file
            * --list     Sample list, format 'SampleID  fqPath' (If PE reads, use comma to separat two fq files)
              --outdir   Output directory, default [./]
Step:
              --steps     Step choosen with comma separating, default [1,2,3,4,5,6,7,8,9,10]
                          1. Inspect raw data
                          2. Filter adapter contained, high Ns, and low-quality sequences
                          3. Align clean reads to the reference genome
                          4. Quantify expression level based on ReadCounts/FPKM/TPM
                          5. Merge ReadCounts/FPKM/TPM of samples to matrixes
                          6. Clustering cells by Seurat
                          7. Detection of alternative splicing events
                          8. SNP calling
                          9. Detection of fusion genes
                          10. Detection of RNA-editing events
Func:
              --template  Generate default configure file
              --help      Print this help information
```
## Example
To run scICW for single-cell RNA sequencing data analysis, the config.txt should be generated at the first place and modified to your needs.
```perl
perl scICW.v1.pl --template
```
Given a sample list ``` sample.list```, with sample id and /path/to/sample.fq.gz a line, you can run scICW with the following command:
```perl
perl scICW.v1.pl --conf config.txt --list sample.list --outdir ./
```