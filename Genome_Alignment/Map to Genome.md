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

    minimap2 -ax map-hifi -t 50 /public/groups/meyerlab/eseal/refgenome/refgen_files/GCF_029215605.1_mMirAng1.0.hap1_genomic.fna \ /public/groups/meyerlab/eseal/refgenome/HiFI_data/SRR25478317.fastq | \ samtools sort -@50 -m 8G -O BAM -o SRR25478317_eseal_sorted.bam - samtools index SRR25478317_eseal_sorted.bam

Because I have high-quality PacBio HiFi data and a map rate of 99% (after flagstat), I don't need to worry about trimming for unmapped reads. HURRAY! Now, we have our sorted BAM file ready for any additional downstream analysis. :) 

Coverage = 35X. PERFECT!!!!!!!!! SHOULD HAVE DONE THIS IS THE FIRST PLACE. Thank you reviewer 2!!!!

03/20/2025: So for some reason the scaffold names (SCAF_...) didn't transfer to the BAM and thus VCF file created. Here is what I did to remedy that: 

1.) Created a tab-deliminated file with the chromosome names in the bam file on one side and the new scaffold names on the right.  
2.) Used samtools view to create a sam file with the old header names 
    
    samtools view -H SRR25478317_eseal_sorted.bam > old_header.sam
3.) Then, used samtools reheader to generate a new bam file with the proper scaffold names

    samtools reheader old_header.sam SRR25478317_eseal.bam > SRR25478317_eseal_sorted_rn.bam
