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

  
