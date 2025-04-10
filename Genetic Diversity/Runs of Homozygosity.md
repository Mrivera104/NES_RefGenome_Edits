# Calculate Runs of Homozygosity
I will be calculating runs of homozygosity using bcftools, just like I did before. 

    bcftools roh -G30 --AF-dflt 0.05 /public/groups/meyerlab/eseal/refgenome/variant_call/SRR25478317_eseal_renamed.vcf.gz.vcf -o SRR25478317_eseal_roh

Next, I isolate only the information I want from the generated file - that is, sample name, chromosome, start position of the ROH, end position, and length (in bp):

    grep "RG" SRR25478317_eseal_roh | cut -f 2,3,4,5,6 > SRR25478317_eseal_roh.txt

Here is the R code I used to visualize results and calculate SROH, NROH, LROH, and FROH

```{r}
# Load the required packages
library(tidyverse)
library(ggrepel)
library(readr)
library(ggplot2)
library(dplyr)


# Set the working directory to where the data files are located
setwd("C:/Users/Millie/Desktop/eseal_popgen") #Change the path to your own folder

# Read the ROH data file
roh_data <- read_delim(file = "SRR25478317_eseal_roh.txt", 
                          delim = "\t", 
                          skip = 1,
                          col_names = c("sample", "scaf", "start", "end", "length (bp)"))

# Subset the data to only include SCAF_1 to SCAF_17 and Length (bp) >= 100,000
roh_data_subset <- roh_data %>% 
  filter(scaf %in% paste0("SCAF_", 1:17)) %>%
  filter(`length (bp)` >= 100000)

roh_mid <- roh_data %>%
  filter(scaf %in% paste0("SCAF_", 1:17)) %>%
  filter(`length (bp)` >= 100000 & `length (bp)` <= 1000000)

# Plot ROH length distribution
ROH_displot <- ggplot(roh_data_subset, aes(x = `length (bp)`)) +
  geom_histogram(binwidth = 10000, fill = "dodgerblue") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  labs(x = "Length (bp)", y = "Count")

ggsave(filename = paste0('Mirounga ROH Distribution Plot'), plot = ROH_displot, 
       device = 'png', dpi = 600, units = c('cm'), width = 28, height = 18, 
       path = "C:/Users/Millie/Desktop/", bg = "white")

# Assuming you already know the total length of the genome
total_genome_length <- 2430321998 

# Compute NROH and SROH
eseal_nroh <- roh_data_subset %>% 
  summarize(NROH = n())

eseal_sroh <- roh_data_subset %>% 
  summarize(SROH = sum(`length (bp)`))

# Compute FROH
eseal_froh <- eseal_nroh %>% 
  bind_cols(eseal_sroh) %>% 
  mutate(FROH = (SROH / total_genome_length) * 100)

# Create a table with NROH, SROH, and FROH for each sample
summary_table <- eseal_froh

# Display the table
print(summary_table)

# sum of all roh greater than 1Mbp 
LROH <- sum(roh_RG_data_subset$`length (bp)`[roh_data_subset$`length (bp)` > 1000000]) 
Mbp_percent <- (LROH/total_genome_length) * 100
mean_SROH <- mean(roh_data_subset$`Length (bp)`)

# percent of genome between 50kbp and 1Mbp
sum_mid <- sum(roh_mid$`length (bp)`)
mid_percent <- (sum_mid/total_genome_length) * 100

# percent of genome in LROH 
print((LROH/total_genome_length) * 100)
```
