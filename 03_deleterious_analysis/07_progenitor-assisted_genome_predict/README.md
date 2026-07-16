## predicted progenitor-assisted genome
Two-script pipeline for analyzing deleterious mutation burden in 100-kb windows:
1.  **DMV6_100K.SNP.window.deleteriousBurden-Min_RMBothMajor-RmSelectiveRegion.R:** Identifies genomic regions where the *minimum* deleterious burden across *Candolleanum* (CND) accessions is lower than in inbreds (A626/E463).
2.  **DMV6_100K.SNP.window-Min_RMBothMajor-RmSelectiveRegion-Min20Founder.R:** Implements a greedy selection algorithm to identify the most influential CND samples ($k=2$–$51$) by maximizing the reduction of the aggregate burden.
