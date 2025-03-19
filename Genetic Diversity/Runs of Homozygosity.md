# Calculate Runs of Homozygosity
I will be calculating runs of homozygosity using bcftools, just like I did before. 

    bcftools roh -G30 --AF-dflt 0.1 /public/groups/meyerlab/eseal/refgenome/variant_call/SRR25478317_eseal_varfilter.vcf -o SRR25478317_eseal_roh
