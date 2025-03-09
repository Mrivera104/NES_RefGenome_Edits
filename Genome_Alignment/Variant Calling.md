# Variant Calling
I am going to be using GATK to perform variant calling on the mapped HiFi data. My BAM is sorted and indexed and is all ready to go. 

Before HaplotypeCaller, I need to process my BAM file by running MarkDuplicates and AddOrReplaceReadGroups. Here is some bash script to help me do so: 

```console
#!/bin/bash

# Define input and output directories
INPUT_DIR=/public/groups/meyerlab/eseal/refgenome 
OUTPUT_DIR="${INPUT_DIR}/variant_call"

# Define the input BAM file
BAM_FILE="/public/groups/meyerlab/eseal/refgenome/your_bam_file.bam"
BASENAME=$(basename "$BAM_FILE" .bam)
MARKDUP_BAM="${OUTPUT_DIR}/${BASENAME}_dupMarked.bam"
METRICS_FILE="${OUTPUT_DIR}/${BASENAME}_dupMetrics.txt"
READGROUP_BAM="${OUTPUT_DIR}/${BASENAME}_dupMarked_RG.bam"

# Run GATK MarkDuplicates
echo "Running MarkDuplicates on $BAM_FILE..."
gatk MarkDuplicates \
  -I "$BAM_FILE" \
  -O "$MARKDUP_BAM" \
  -M "$METRICS_FILE"

# Run GATK AddOrReplaceReadGroups
echo "Running AddOrReplaceReadGroups on $MARKDUP_BAM..."
gatk AddOrReplaceReadGroups \
  -I "$MARKDUP_BAM" \
  -O "$READGROUP_BAM" \
  -ID 1 \
  -LB lib1 \
  -PL PACBIO \
  -PU unit1 \
  -SM SRR25478317

```

After all this, we should have BAM files that are ready to be processed. To perform variant calling, we will first use haplotypecaller to create a genotype VCF file, call genotypes and filter to get a proper VCF file we can use for downstream analyses. 
```console
#!/bin/bash

# Define input and output directories
INPUT_DIR="/public/groups/meyerlab/eseal/refgenome/variant_call"
OUTPUT_DIR="${INPUT_DIR}"

# Define input files
REFGENOME="GCF_029215605.1_mMirAng1.0.hap1_genomic.fna"
BAM_FILE="SRR25478317_eseal_sorted_dupMarked_RG.bam"
BASENAME=$(basename "SRR25478317_eseal")
GVCF_FILE="${OUTPUT_DIR}/${BASENAME}.g.vcf.gz"
VCF_FILE="${OUTPUT_DIR}/${BASENAME}.vcf.gz"

# Run GATK HaplotypeCaller
echo "Running HaplotypeCaller on $BAM_FILE..."
gatk HaplotypeCaller \
  -I "$BAM_FILE" \
  -R "$REFGENOME" \
  -ERC GVCF \
  -O "$GVCF_FILE"

if [ $? -eq 0 ]; then
  echo "Finished HaplotypeCaller for $BAM_FILE."
else
  echo "Error running HaplotypeCaller on $BAM_FILE."
  exit 1
fi

# Index VCF file
tabix -p vcf "$GVCF_FILE"

# Run GATK GenotypeGVCFs
gatk GenotypeGVCFs \
  -R "$REFGENOME" \
  -V "$GVCF_FILE" \
  -O "$VCF_FILE"

if [ $? -eq 0 ]; then
  echo "Finished GenotypeGVCFs for $GVCF_FILE."
else
  echo "Error running GenotypeGVCFs for $GVCF_FILE."
  exit 1
fi

# Run GATK VariantFiltration
gatk VariantFiltration \
  -R "$REFGENOME" \
  -V "$VCF_FILE" \
  -O "${OUTPUT_DIR}/${BASENAME}_varfilter.vcf" \
  --filter-name "Low_depth10" \
  --filter-expression "DP < 10"

if [ $? -eq 0 ]; then
  echo "Finished VariantFiltration for $VCF_FILE."
else
  echo "Error running VariantFiltration for $VCF_FILE."
  exit 1
fi
```
