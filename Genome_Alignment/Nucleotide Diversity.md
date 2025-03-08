# Calculate Genome-Wide Heterozygosity

I'm going to go about this in two ways. The first is calculating genome-wide heterozygosity like I did in my first go-around. The second is to calculate using ANGSD. Maybe this time, it will be different. 

I won't subset by the largest scaffolds this time; I want to include them all.

# Calculate Genome-Wide Heterozygosity using ANGSD

We can use ANGSD to calculate nucleotide diversity. Here is the explanation: "The heterozygosity is the proportion of heterozygous genotypes. This is in some sense encapsulated in the theta estimates." http://www.popgen.dk/angsd/index.php/Heterozygosity

From ChatGPT: "The heterozygosity value provided by ANGSD can be interpreted as nucleotide diversity, reflecting the average genetic variation at the nucleotide level across the genome."

    angsd -P 8 -i SRR25478317_eseal_sorted.bam -anc /public/groups/meyerlab/eseal/refgenome/refgen_files/GCF_029215605.1_mMirAng1.0.hap1_genomic.fna -dosaf 1 -gl 1 -out SRR25478317_eseal_angsdput -C 50 -ref /public/groups/meyerlab/eseal/refgenome/refgen_files/GCF_029215605.1_mMirAng1.0.hap1_genomic.fna -minQ 20 -minmapq 30
    realSFS -fold 1 -P 80 SRR25478315_eseal_mapped_angsdput.saf.idx > SRR25478315_eseal_mapped_est.ml



