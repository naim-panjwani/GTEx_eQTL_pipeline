# GTEx_eQTL_pipeline

This repository analyzes eQTLs for brain and nerve tissues per sex using the [GTEx eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) docker image.  


# Subset samples and tissues desired

`00-setup.sh` should be run first to setup folders, symbolic links and cleaner covariate files, PLINK files, etc.

`phenotypes_eda.ipynb` shows the reasoning on how to list out the tissue sample ID's per brain and tibial nerve tissues and per sex

`01-run_eqtl_pipeline.sh` goes through the [GTEx eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) as closely as possible for the sample subsets that were generated in `phenotypes_eda.ipynb` step

`01-submit_eqtl_pipeline` is to submit each of the sample subsets generated as a job on the hpc

eQTL analyses were adjusted for:  

  - Top 5 genotyping principal components.  
  - A set of covariates identified using the Probabilistic Estimation of Expression Residuals (PEER) method (Stegle et al., PLoS Comp. Biol., 2010 ), calculated for the normalized expression matrices (described below). For eQTL analyses, the number of PEER factors was determined as function of sample size (N): 15 factors for N<150, 30 factors for 150≤ N<250, 45 factors for 250≤ N<350, and 60 factors for N≥350, as a result of optimizing for the number of eGenes discovered. For sQTL analyses, 15 PEER factors were computed for each tissue.  
  - Sequencing platform (Illumina HiSeq 2000 or HiSeq X).  
  - Sequencing protocol (PCR-based or PCR-free).  
  - Sex. Although this was dropped by `combined_covariates.py` due to collinearity (since each subset is one sex).  
  - Age (note: this is adjusted for in the public analyses).  

