# README: Genome assembly and annotation of a deep-diving pinniped, the northern elephant seal (*Mirounga angustirostris*)
Genomic Analyses for the NES Reference Genome Manuscript 

## Overview

This repository contains scripts used for genomic analyses using the reference genome of the northern elephant seal (*Mirounga angustirostris*) created through the California Conservation Genomics Project (CCGP). Analyses conducted in this paper include repeat masking, genome alignment of PacBio Hifi long read data, variant calling, estimation of genome-wide heterozygosity, and variant annotation. Scripts used in the de novo assembly of the NES reference genome are provided by the CCGP and included in the manuscript. 

The repository is organized to maximize reproducibility. All scripts use relative paths where possible. Input data required to reproduce R Studio plots are included.

# Description of the data

 We use the raw FASTA accompanying annotation file downloaded through NCBI BioProject accession number PRJNA937321. GenBank assembly for the raw FASTA file is GenBank assembly GCA_029215605.1. The annotation file used was GCF_029215605.1_mMirAng1.0.hap1_genomic.gff.gz downloaded through NCBI FTP.

### How to download the data

The following code was used to download the raw FASTA file and PacBio HiFi long read data from NCBI. 

```
datasets download genome accession GCF_02921506.1 --include genome
```
```
prefetch --max-size 100G SRR25478317
fasterq-dump --threads 8 SRR25478317.sra
```
