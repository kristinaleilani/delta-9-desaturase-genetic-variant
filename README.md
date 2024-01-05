# delta-9-desaturase-genetic-variant
Code associated with the paper: A genetic variant of delta-9 desaturase is associated with cold hardiness in a Great Barrier Reef coral

Pipelines for processing raw sequencing data (for 2bRAD, TagSeq, and whole-genome phased datasets) are in Processing sequences.sh. All processing was done on Lonestar 6  at the Texas Advanced Computing Center.

The following R scripts reproduce figures from the manuscript:

Figure 1:
Admixture.R = Admixture barplot
Bayescan_plots.R = Fst Manhattan plot 
Sambada.R = Map of allele turnover across the GBR

Figure 2:
Gene_expression.R = DESeq2 volcano plots

Figure 3:
Gene_expression.R = Comparison of delta-9 desaturase expression between Orpheus and Keppel corals in a reciprocal transplant experiment.
Growth.R = Comparison of growth between Orpheus and Keppel corals in a reciprocal transplant experiment

Figure 4:
chromomap.R = Map of delta-9 desaturase genes on chromosome 9, plus significant loci



For questions/clarification, please to reach out to Kristina Black (kblack@utexas.edu).

