# Calculating the amount of LoF variants and how many LoF variants there are per chromosome

For creating a loss of function variant rate per chromosome figure in R, I need to create a TSV file with all the necessary information. For that, I used this bash script: 

```
SnpSift extractFields SRR25478317_eseal_output_homsites_subset.ann.vcf \
"CHROM" "POS" "REF" "ALT" \
"ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \
"ANN[*].FEATURE_ID" "ANN[*].BIOTYPE" | \
grep -E 'stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_lost|exon_loss_variant' > LOF_variants.tsv
```
This matches what I did for the g:Profiler analysis. 

My R script for turning this into a LoF variant rate figure is as follows: 

```
setwd("C:/Users/Millie/Desktop/eseal_refgenome/manuscript_edits")

library(readr)
library(ggplot2)
library(dplyr)
library(viridis)

# Load and clean data
raw_data <- readr::read_tsv("LOF_variants2.tsv", col_names = FALSE)
colnames(raw_data) <- c("chrom", "position", "ref", "alt")

# Chromosome name mapping
chrom_order <- c("NW_026991709.1", "NW_026991710.1", "NW_026991711.1", "NW_026991712.1", 
                 "NW_026991713.1", "NW_026991714.1", "NW_026991715.1", "NW_026991716.1",
                 "NW_026991717.1", "NW_026991718.1", "NW_026991719.1", "NW_026991720.1",
                 "NW_026991721.1", "NW_026991722.1", "NW_026991723.1", "NW_026991724.1",
                 "NW_026991725.1")

# Plot 1: LoF count per chromosome (barplot)
variants_per_chrom <- as.data.frame(table(raw_data$chrom)) 
colnames(variants_per_chrom) <- c("chrom", "count")

variants_per_chrom <- variants_per_chrom %>%
  mutate(chrom = recode(chrom, 
                        NW_026991709.1 = 'Scaf_1', NW_026991710.1 = 'Scaf_2', 
                        NW_026991711.1 = 'Scaf_3', NW_026991712.1 = 'Scaf_4', 
                        NW_026991713.1 = 'Scaf_5', NW_026991714.1 = 'Scaf_6', 
                        NW_026991715.1 = 'Scaf_7', NW_026991716.1 = 'Scaf_8', 
                        NW_026991717.1 = 'Scaf_9', NW_026991718.1 = 'Scaf_10', 
                        NW_026991719.1 = 'Scaf_11', NW_026991720.1 = 'Scaf_12', 
                        NW_026991721.1 = 'Scaf_13', NW_026991722.1 = 'Scaf_14', 
                        NW_026991723.1 = 'Scaf_15', NW_026991724.1 = 'Scaf_16', 
                        NW_026991725.1 = 'Scaf_17'))

chrom_length = c(215935216, 
           205257464, 
           195234058,
           190123178,
           180400003,
           154989474,
           154172358,
           146104457,
           144560927,
           138646194,
           127380783,
           109230848,
           108890317,
           95463112,
           93938456,
           60166985,
           57920945)

variants_per_chrom$length <- chrom_length
variants_per_chrom$variant_length <- variants_per_chrom$length / variants_per_chrom$count
variants_per_chrom$variants_per_mb <- 1e6 / variants_per_chrom$variant_length

LoF_var_rate_plot <- ggplot(variants_per_chrom, aes(x = chrom)) +
  geom_bar(aes(y = count), stat = "identity", fill = "#8dccc5", width = 0.6) +
  geom_line(aes(y = variants_per_mb * (max(count) / max(variants_per_mb)), group = 1), 
            color = "tomato2", size = 0.7) +
  geom_point(aes(y = variants_per_mb * (max(count) / max(variants_per_mb))), 
             color = "tomato2", size = 2) +
  scale_y_continuous(
    name = "Number of Variants",
    sec.axis = sec_axis(~ . / (max(variants_per_chrom$count) / max(variants_per_chrom$variants_per_mb)), name = "Variants per Mb"),
    expand = c(0, 0), limits = c(0, 500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x=element_blank(),
        axis.title.y= element_text(size = 16, color = "black", family = "Arial"),
        axis.text.y.right = element_text(size = 10, margin = margin(15, 15, 15, 15)),
        axis.text.y.left = element_text(size = 10, margin = margin(15, 15, 15, 15)))

print(LoF_var_rate_plot)

ggsave(filename = 'LoF_variant_rate.png', plot = LoF_var_rate_plot, 
       device = 'png', dpi = 600, units = c('cm'), width = 28, height = 18, 
       path = "C:/Users/Millie/Desktop/", bg = "white")

print(mean(variants_per_chrom$count))
print(mean(variants_per_chrom$variants_per_mb))
```

