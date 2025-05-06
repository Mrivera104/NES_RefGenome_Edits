# Inferring historical demography using MSMC2

I was advised to use MSMC2 instead of PSMC, so that's what I'm going to (try to) do. 

First, we need to generate the required mask files: one mappability bed file, and one bed file containing call information. 

from the MSMC github: 

CALL INFORMATION BED FILE:
- MSMC is a Hidden Markov Model which uses the density of heterozygous sites (1/0 genotypes) to estimate the time to the most recent common ancestor. However, for a density you need not only a numerator but also a denominator, which in this case is the number of non-heterozygous sites, so typically homozygous reference alleles. Those are not part of this VCF file, for efficiency reasons. That's why we have a Mask-file for each sample, that gives information in which regions in the genome could be called.
- Regions with not enough coverage or too low quality will be excluded.

Looking at this now, since I have all sites emitted, I might not actually need this mask. But whatevs, I will make it anyways and run it once with and once without it. 

 MAPPABILITY BED FILE: 
- This mask defines regions in the reference genome in which we trust the mapping to be of high quality because the reference sequence is unique in that area.

# Step 1: Generate Mask Files
I am going to use samtools depth to create a bed file with callable sites. 

    samtools depth -a /public/groups/meyerlab/eseal/refgenome/SRR25478317_eseal_sorted_subset.bam | \
      awk '$3 >= 10 {print $1"\t"$2-1"\t"$2}' > callable_mask.bed

For the mappability bed file, I am using a software called GenMap. 

    genmap index -F GCF_029215605.1_mMirAng1.0.hap1_genomic.fna -I genmap_index_dir
    genmap map -K 100 -E 2 -I genmap_index_dir -O . -t 4
- `K 100`: read length
- `E 2`: allow 2 mismatches
- `I genmap_index_dir`: write to a separate directory
- `O`: output directory
- `-t`: how many cores to use

Result is a .wig file showing mappability scores from 0 to 1.
