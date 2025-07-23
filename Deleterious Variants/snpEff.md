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

# Step 4: Filter Variants from Annotated VCF 

For this first analysis, I'm trying to get a unique gene name list of high impact variants. These are the effects I've chosen as high impact, putative loss of function effects: 

- stop gained
- frameshift variant
- splice acceptor variant
- splice donor variants
- start lost
- stop lost
- exon loss variant

Here is code I used to generate this gene list: 

```
SnpSift filter \
  "(ANN[*].IMPACT = 'HIGH') & ((ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'start_lost') | (ANN[*].EFFECT has 'stop_lost') | (ANN[*].EFFECT has 'exon_loss_variant'))" \
  SRR25478317_eseal_output_homsites_subset.ann.vcf | \
SnpSift extractFields - -s "\t" "ANN[*].GENE" | \
grep -v "^#" | sort | uniq > SRR25478317_LoF_HighImpact_genes_2.txt
```

For creating a loss of function variant rate per chromosome figure in R, I need to create a TSV file with all the necessary information. For that, I used this bash script: 

```
SnpSift extractFields SRR25478317_eseal_output_homsites_subset.ann.vcf \
"CHROM" "POS" "REF" "ALT" \
"ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \
"ANN[*].FEATURE_ID" "ANN[*].BIOTYPE" | \
grep -E 'stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_lost|exon_loss_variant' > LOF_variants.tsv
```




