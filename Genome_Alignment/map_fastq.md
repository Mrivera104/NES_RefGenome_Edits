# Mapping HiFi long-read data to CCGP reference genome
** 03/06/2025 **
This is the beginning of the edits I'm doing for the reference genome paper. I was advised to use the higher quality long-read data instead of the Omni-C short-read data for my downstream analyses. So I'm going to do that! 

** All of this is done in a conda environment making it easy to download required packages **

I downloaded the data I needed from NCBI.


> Downloading reference genome (hap1)

    datasets download genome accession GCF_02921506.1 --include genome

> Download HiFi long-read data

    prefetch --max-size 100G SRR25478317
    fasterq-dump --threads 8 SRR25478317.sra

Now we are pretty much ready to map the long-read fastq files to the reference genome. Unlike before, I have to use Minimap2 for alignment purposes, since BWA-MEM works only for short-read data. 

    minimap2 -a -t 50 /public/groups/meyerlab/eseal/refgenome/refgen_files/GCF_029215605.1_mMirAng1.0.hap1_genomic.fna /public/groups/meyerlab/eseal/refgenome/HiFI_data/SRR25478317.fastq | samtools sort -@50 -O BAM -o SRR25478317_eseal_sorted.bam -


