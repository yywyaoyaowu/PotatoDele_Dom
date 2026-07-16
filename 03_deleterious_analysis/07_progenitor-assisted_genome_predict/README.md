## Deleterious Mutation Burden Analysis

Two core R scripts for quantifying and dissecting genetic load across the potato genome.

### 1. Genome-wide Burden Matrix Construction (`generate_burden_matrix.R`)
Constructs a genome-wide deleterious burden matrix from VCF frequency data and conservation-based genotype effects.
*   **Sliding Windows:** Partitions the genome into 100-kb windows (excluding selective sweeps).
*   **Burden Calculation:** Calculates weighted deletion burdens (GERP scores × genotype) for ~500 accessions.
*   **Baseline Output:** Generates `DMRef_100KSNP.window_AllDeleBurden...bed`, the master table for downstream analyses.

### 2. Comparative Burden & Core Accession Identification (`analyze_burden.R`)
Identifies genomic regions and specific accessions driving genetic load in *Candolleanum* (CND).
*   **Elite Comparison:** Filters windows where the **minimum** burden across CND accessions is lower than in elite cultivars (A626 and E463).
*   **Greedy Optimization:** Implements a greedy selection algorithm ($k=2$–$51$) to identify the minimal subset of CND samples that maximizes the reduction of aggregate burden.
*   **Output:** Produces BED files of candidate regions and CSV/TXT reports detailing the contribution of specific accessions to the overall genetic load.
