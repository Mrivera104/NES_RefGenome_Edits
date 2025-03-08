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

if [ $? -eq 0 ]; then
  echo "Finished MarkDuplicates for $BAM_FILE. Output: $MARKDUP_BAM"
else
  echo "Error running MarkDuplicates on $BAM_FILE."
  exit 1
fi

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

if [ $? -eq 0 ]; then
  echo "Finished AddOrReplaceReadGroups for $MARKDUP_BAM. Output: $READGROUP_BAM"
else
  echo "Error running AddOrReplaceReadGroups on $MARKDUP_BAM."
  exit 1
fi

echo "Processing completed successfully. Output files are in $OUTPUT_DIR."
```
