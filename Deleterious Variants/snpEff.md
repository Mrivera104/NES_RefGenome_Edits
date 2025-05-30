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


    java -jar SnpSift.jar filter "ANN[0].EFFECT has 'missense_variant' | ANN[0].EFFECT has 'frameshift_variant'" SRR25478317_eseal_output_homsites_subset.ann.vcf > deleterious_variants.vcf


    java -jar SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" SRR25478317_eseal_output_homsites_subset.ann.vcf > high_impact_output.vcf


    java -jar SnpSift.jar filter "(QUAL > 30) & (DP > 10)" high_impact_output.vcf > high_impact_qual_filter_output.vcf

