# README: Genome assembly and annotation of a deep-diving pinniped, the northern elephant seal (*Mirounga angustirostris*)
Genomic Analyses for the NES Reference Genome Manuscript 

## Overview of the data

This repository contains scripts used for genomic analyses using the reference genome of the northern elephant seal (*Mirounga angustirostris*) created through the California Conservation Genomics Project (CCGP). Analyses conducted in this paper include repeat masking, genome alignment of PacBio Hifi long-read data, variant calling, estimation of genome-wide heterozygosity, and variant annotation. Scripts used in the de novo assembly of the NES reference genome are provided by the CCGP. The github repository to all genome assembly scripts included in the published manuscript. 

The repository is organized to maximize reproducibility. All scripts use relative paths where possible. Input data required to reproduce R Studio plots are included.

### Data sources
- Reference genome: GenBank assembly GCA_029215605.1 (BioProject PRJNA937321)
- PacBio HiFi long-reads: NCBI SRA accession SRR25478317 (BioProject PRJNA998853)
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

## Alignment of long-read data to genome
- Tools: Minimap2 v2.29, Samtools v1.11
- Purpose: Align long-reads to reference genome, sort and index BAM files.

## Variant Calling
- Tool: GATK v3.8
- Purpose: Detect SNPs and indels. Pipeline includes duplicate marking, read group assignment, and variant filtering.

## Estimation of genome-wide heterozygosity
- Tool: ANGSD v0.94
- Purpose: Calculate heterozygosity across the genome.

## Historical demography
- Tools: BCFtools v1.14, PSMC v0.6.5
- Purpose: Infer effective population size changes over time.

## Variant annotation
- Tools: SnpEff v5.2, SnpSift v5.2
- Purpose: Annotate variants, classify loss-of-function (LoF) mutations.

# Contents of this repository

## Bash scripts
- repeat_mask.sh - runs RepeatMasker on the genome assembly to identify and annotate repetitive elements.
- genome_alignment.sh - aligns sequencing reads to the reference genome using appropriate mappers.
- variant_calling.sh – pipeline for read alignment and variant calling.
- angsd_heterozygosity.sh - estimates genome-wide heterozygosity using ANGSD.
- historical_demography.sh - generates input files and runs PSMC to infer historical population size changes.
- variant_annotation.sh - annotates genetic variants with predicted functional effects using standard annotation tools.
- extract_lof_variants.sh - filters annotated variants to extract those predicted to cause loss-of-function.
- gene_name_extraction - pulls gene names associated with LoF variants for downstream analysis.

## R scripts
- nucleotide_diversity_comparison.R – marine mammal genome-wide heterozygosity comparison plot.
- lof_variant_plotting.R – visualization of LoF variant density and rate across scaffolds.

## Supporting data
- nuc_diversity_marinemammals.csv – compiled heterozygosity values across marine mammals.
- LOF_variants.tsv – annotated LoF variants in the elephant seal genome.
- LoF_gene_list.txt - list of gene names associated with LoF variants. 

# Citations
- Bao, W., Kojima, K. K., & Kohany, O. (2015). Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA, 6, 1-6.
- Cingolani, P. (2012). Variant annotation and functional prediction: SnpEff. In Variant calling: methods and protocols (pp. 289-314). New York, NY: Springer US.
- Cingolani, P., Patel, V. M., Coon, M., Nguyen, T., Land, S. J., Ruden, D. M., & Lu, X. (2012). Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift. Frontiers in genetics, 3, 35.
- Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10(2), giab008.
- Fiedler, P. L., Erickson, B., Esgro, M., Gold, M., Hull, J. M., Norris, J. M., ... & Shaffer, H. B. (2022). Seizing the moment: the opportunity and relevance of the California Conservation Genomics Project to state and federal conservation policy. Journal of Heredity, 113(6), 589-596.
- Korneliussen, T. S., Albrechtsen, A., & Nielsen, R. (2014). ANGSD: analysis of next generation sequencing data. BMC bioinformatics, 15, 1-13.
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100.
- Li, H., & Durbin, R. (2011). Inference of human population history from individual whole-genome sequences. Nature, 475(7357), 493-496.
- Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., ... & Banks, E. (2017). Scaling accurate genetic variant discovery to tens of thousands of samples. BioRxiv, 201178.
- Smit, A. F. A., Hubley, R., & Green, P. (2015). RepeatMasker Open-4.0. 2013–2015.
