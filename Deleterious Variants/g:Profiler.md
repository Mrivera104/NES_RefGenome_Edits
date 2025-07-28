# Performing Gene Enrichment Analysis of LoF variants using g:Profiler

We can use the output from SnpEff to plug into g:Profiler. g:Profiler is a public web server for characterising and manipulating gene lists. 

# g:GOST 

From the website: 

"g:GOSt performs functional profiling of gene lists using various kinds of biological evidence. The tool performs statistical enrichment analysis to find over-representation of information from Gene Ontology terms, biological pathways, regulatory DNA elements, human disease gene annotations, and protein-protein interaction networks."
"g:GOSt uses Fisher's one-tailed test, also known as cumulative hypergeometric probability, as the p-value measuring the randomness of the intersection between the query and the ontology term. The p-value represents the probability of the observed intersection plus probabilities of all larger, more extreme intersections."

# Filter High Impact LoF Variants from Annotated VCF 

For this analysis, I'm trying to get a unique gene name list of high impact variants. These are the effects I've chosen as high impact, putative loss of function effects: 

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

Then I just plug into g:Profiler, mapping the gene list to the American mink annotation. I had to manually go in there and make sure exon:ectect thingies were gone. g:Profiler doesn't care if you put in duplicates of a gene name. 
