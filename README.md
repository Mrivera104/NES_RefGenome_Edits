# README: Genome assembly and annotation of a deep-diving pinniped, the northern elephant seal (*Mirounga angustirostris*)
Genomic Analyses for the NES Reference Genome Manuscript 

## Overview

This Dryad repository contains scripts, raw input files, and processed output files used for genomic analyses of the northern elephant seal (*Mirounga angustirostris*). Analyses include 

The repository is organized to maximize reproducibility. All scripts use relative paths where possible, and all input data required to reproduce plots are included.

## Data Types
- Raw sequencing data: PacBio HiFi reads for a single individual.
- Reference genome: Haploid genome assembly with accompanying GFF3 annotations.
- Aligned data: BAM files generated from HiFi reads aligned to the reference genome.
- Variant calls: GVCF, VCF, filtered VCF, and subset VCF files from GATK.
- Plots: Visualizations of LoF variant counts and genome-wide heterozygosity in marine mammals.
- Supporting data: CSV files used for heterozygosity and LoF analyses.
- Scripts: Bash scripts for alignment, variant calling, RepeatMasker, PSMC demography, and R scripts for plotting.
