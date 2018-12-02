#!/bin/bash


tmp=181201

source 00-set_variables.sh
if [ ! -d "$step07_outputdir" ]; then mkdir $step07_outputdir; fi

famfile="${step04_outputdir}/${tissue_of_interest}_WGS_genotypes_bplink.fam" # this is the WGS samples only; no Omni here
peers_file="${step05_outputdir}/${ExpressionFile%.gct.gz}_${tissue_of_interest}_normalized.PEER_covariates.txt"
pca_file="${step06_outputdir}/${tissue_of_interest}_PC-AiR_eigenvectors.txt"
subj_sample_attr_file="${PhenoDir}/Covariates.txt"

Rscript 07-Generate_covariates.R $famfile \
				 $subj_sample_attr_file \
				 $pca_file \
				 $peers_file \
				 ${step07_outputdir}/${tissue_of_interest}_covariates.txt

