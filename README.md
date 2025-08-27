# README: Genome assembly and annotation of a deep-diving pinniped, the northern elephant seal (*Mirounga angustirostris*)
Genomic Analyses for the NES Reference Genome Manuscript 

## Overview of the data

This repository contains scripts used for genomic analyses using the reference genome of the northern elephant seal (*Mirounga angustirostris*) created through the California Conservation Genomics Project (CCGP). Analyses conducted in this paper include repeat masking, genome alignment of PacBio Hifi long read data, variant calling, estimation of genome-wide heterozygosity, and variant annotation. Scripts used in the de novo assembly of the NES reference genome are provided by the CCGP and included in the manuscript. 

The repository is organized to maximize reproducibility. All scripts use relative paths where possible. Input data required to reproduce R Studio plots are included.

### Data sources

- Reference genome: GenBank assembly GCA_029215605.1 (BioProject PRJNA937321)
- PacBio HiFi long reads: NCBI SRA accession SRR25478317 (BioProject PRJNA998853)
- Genome annotation: GCF_029215605.1_mMirAng1.0.hap1_genomic.gff.gz (NCBI FTP)

### Download commands

```
# Reference genome
datasets download genome accession GCF_02921506.1 --include genome
```
```
# PacBio HiFi reads
prefetch --max-size 100G SRR25478317
fasterq-dump --threads 8 SRR25478317.sra
```
```
# Annotation file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/215/605/GCF_029215605.1_mMirAng1.0.hap1/GCF_029215605.1_mMirAng1.0.hap1_genomic.gff.gz
```

# Analyses

## Repeat Masking

- Tool: RepeatMasker v4.1.5
- Purpose: Mask interspersed repeats and low-complexity regions in the reference genome.

## Alignment of long read data to genome

- Tools: Minimap2 v2.29, Samtools v1.11
- Purpose: Align long reads to reference genome, sort and index BAM files.

## Variant Calling

- Tool: GATK v3.8
- Purpose: Detect SNPs and indels. Pipeline includes duplicate marking, read group assignment, and variant filtering.

## Estimation of genome-wide heterozygosity

- Tool: ANGSD v0.94
- Purpose: Calculate heterozygosity from sequence data across the genome.

## Historical demography

- Tools: BCFtools v1.14, PSMC v0.6.5
- Purpose: Infer effective population size changes over time.

## Variant annotation

- Tools: SnpEff v5.2, SnpSift v5.2
- Purpose: Annotate variants, classify loss-of-function (LoF) mutations.

# Contents of this repository

## Bash scripts

- gatk_variant_calling_pipeline.sh – pipeline for read alignment and variant calling

## R scripts

- lof_variant_plot.R – visualization of LoF variant density across scaffolds
- nuc_diversity_barplot.R – marine mammal heterozygosity comparison plot

## Supporting data

- nuc_diversity_marinemammals.csv – compiled heterozygosity values across marine mammals
- LOF_variants.tsv – annotated LoF variants in the elephant seal genome

# Sharing/Access information

