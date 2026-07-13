## Recombination Rate Inference (ReLERNN)
We employed ReLERNN (https://github.com/kr-colab/ReLERNN), a deep learning approach utilizing Convolutional Neural Networks (CNNs), to infer fine-scale recombination rates across the genome. 
The pipeline consisted of four stages:
1. **Simulation:** Generated training, validation, and testing datasets using `ReLERNN_SIMULATE`, parameterized with an assumed mutation rate ($\mu=1\times10^{-8}$) and an upper rho-theta ratio ($R=10$). Simulations used 15-kb windows and ran with a fixed random seed (42) for reproducibility.
2. **Training:** Trained the CNN model using `ReLERNN_TRAIN` on the simulated datasets.
3. **Prediction:** Predicted recombination rates for genomic windows in the empirical VCF data using `ReLERNN_PREDICT`.
All analyses were executed within a dedicated Conda environment (`ML_python3.9`) with CUDA 12.2 support for GPU acceleration.

eg: bash ReLERNN.sh Lan_chr01_output Lan_chr01.snp.recode.vcf Solanum_tuberosumDM_chr01.bed


