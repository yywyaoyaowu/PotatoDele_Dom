## Progenitor-Assisted Genome Prediction
Constructs a genome-wide deleterious burden matrix by integrating GERP conservation scores with genotype data across 100-kb SNP windows (whole genome：01_DM_100K.SNP.window.deleteriousBurden.R；Without selective sweeps ：02_DM_100K.SNP.window.deleteriousBurden_WithoutSelectiveSweep.R ). The script calculates the segments and materials with lower burden than inbreds (A626 and E463).
03_donor.R: get the top five donors with least deleterious burden.
