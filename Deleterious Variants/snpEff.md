# Using SnpEff to annotate snps and look for deleterious variants

SnpEff is a tool that annotates and predicts the effects of genetic variants (such as amino acid changes). By creating a custom database using annotation files, the tool will create an annotated version of the input vcf file to have gene effect information. 

SnpEff and SnpSift documentation: http://pcingola.github.io/SnpEff/

# Step 1: Build a custom SnpEff database

Because I am working with a nonmodel organism, I am required to build my own SnpEff database. First, we make a custom directory in the SnpEff data directory. 

    mkdir data/GCF_029215605.1

Then, I copied the reference genome (FASTA) and annotation (GFF3, GTF, or Gencode) files into the new directory.

    wget (FTP link to GFF and GTF files on NCBI) 
    cp /path/to/reference_genome.fasta .

Next, I modified the snpEff.config file to include my organism and its genome. I did this using nano on the config file; I forgot to record the exact thing I did for this but the steps are on the SnpEff website where they tell you how to do this and where in the config file to add your information. Not hard to look this back up if I do it again. 

# Step 2: Build the SnpEff database 

    java -jar snpEff.jar build -gff3 -v GCF_029215605.1

# Step 3: Annotating Variants with SnpEff

    java -Xmx4g -jar snpEff.jar GCF_029215605.1 /public/groups/meyerlab/eseal/refgenome/variant_call/ncbi/SRR25478317_eseal_output_homsites_subset.vcf.gz > SRR25478317_eseal_output_homsites_subset.ann.vcf

# Step 4: Reporting variants in protein-coding sequences, CDS heterozygosity and highest impact per variant

So now I need to figure out just how much of the genome is in protein-coding sequences. This will allow me to calculate the heterozgosity in these CDS, which will obviously be even lower than genome-wide heterozygosity. I can use the gff file for this. 

```
# Define your list of 17 largest scaffolds
scaffolds="NW_026991709.1 NW_026991710.1 NW_026991711.1 NW_026991712.1 NW_026991713.1 NW_026991714.1 NW_026991715.1 NW_026991716.1 NW_026991717.1 NW_026991718.1 NW_026991719.1 NW_026991720.1 NW_026991721.1 NW_026991722.1 NW_026991723.1 NW_026991724.1 NW_026991725.1"

# Step 1: Extract CDS and convert to BED
awk '$3 == "CDS" {print $1"\t"($4-1)"\t"$5}' genes.gff > cds_all.bed

# Step 2: Filter to only top 17 scaffolds
grep -F -w -f <(echo "$scaffolds" | tr ' ' '\n') cds_all.bed > cds_top17.bed

# Step 3: Merge overlapping regions
bedtools sort -i cds_top17.bed | bedtools merge -i - > cds_merged_top17.bed

# Step 4: Sum total CDS length
awk '{sum += $3 - $2} END {print sum}' cds_merged_top17.bed
```

Next, I had to extract the most impactful consequence per variant, instead of just reporting all possible effects for every variant.

    bcftools view -v snps,indels SRR25478317_eseal_output_homsites_subset.ann.vcf | \
      SnpSift extractFields - "ANN[0].IMPACT" | \
      grep -v '^$' | sort | uniq -c | sort -nr

That leaves me with a total of 1,417,101 effects (which matches the number of variants reported by SnpEff). This is the breakdown: 

- 1401772 MODIFIER 98.91% -> outside of a coding region
-  4970 LOW 0.35%
- 4686 HIGH 0.33%
- 3180 MODERATE 0.22%
- 4970 + 4686 + 3180 = 12836, 0.90% in CDS
- 4970+4686+3180 = 12836 in CDS
- *12836/34767161 = 0.000369 heterozygosity in CDS*







