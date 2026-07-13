## Principal Component Analysis (PCA)
Principal Component Analysis was performed using PLINK on the filtered SNP dataset (MAF ≥ 0.05) to visualize population structure. We extracted the top 502 principal components to capture fine-scale genetic variation. The final eigenvector matrix was formatted with appropriate headers (PC1–PC502) and tab-separated values to facilitate visualization.

eigenvec_plot.R: plot eigenvector matrix

eigenval_plot.R: plot eigenvalue matrix 
