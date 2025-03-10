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

Turns out, after adding readgroups and all of that, you shouldn't perform variant calling with GATK ... Apparently the best way to perform variant calling on long-read data is by using DeepVariant. 

"DeepVariant is a deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file." https://github.com/google/deepvariant

BUT DEEPVARIANT DIDN'T WORK FOR ME!!!!!!!!!! TOO COMPLEX!!!! NEEDS ULTRA COMPUTING POWER I DON'T HAVE! Whatever, I updated GATK and tried again. 

```console
#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Define input and output directories
INPUT_DIR="/public/groups/meyerlab/eseal/refgenome/variant_call"
OUTPUT_DIR="${INPUT_DIR}"

# Define input files
REFGENOME="GCF_029215605.1_mMirAng1.0.hap1_genomic.fna"
BAM_FILE="SRR25478317_eseal_sorted_dupMarked_RG.bam"
BASENAME=$(basename "$BAM_FILE" .bam)  # Extract filename without extension
GVCF_FILE="${OUTPUT_DIR}/${BASENAME}.g.vcf.gz"
VCF_FILE="${OUTPUT_DIR}/${BASENAME}.vcf.gz"
FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_varfilter.vcf"

# Run GATK HaplotypeCaller with 50 threads
echo "Running HaplotypeCaller on $BAM_FILE with 50 threads..."
gatk HaplotypeCaller \
  -I "$BAM_FILE" \
  -R "$REFGENOME" \
  -O "$GVCF_FILE" \
  --emit-ref-confidence GVCF \
  --native-pair-hmm-threads 50

echo "Finished HaplotypeCaller."

# Index GVCF file
echo "Indexing GVCF file..."
tabix -p vcf "$GVCF_FILE"

# Run GATK GenotypeGVCFs
echo "Running GenotypeGVCFs..."
gatk GenotypeGVCFs \
  -R "$REFGENOME" \
  -V "$GVCF_FILE" \
  -O "$VCF_FILE"

echo "Finished GenotypeGVCFs."

# Run GATK VariantFiltration
echo "Running VariantFiltration..."
gatk VariantFiltration \
  -R "$REFGENOME" \
  -V "$VCF_FILE" \
  -O "$FILTERED_VCF" \
  --filter-name "Low_depth10" \
  --filter-expression "DP < 10"

echo "Finished VariantFiltration."
echo "Pipeline completed successfully!"
```
