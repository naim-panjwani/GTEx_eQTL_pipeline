#!/bin/bash

tmp=181201

source 00-set_variables.sh
if [ ! -d "$step06_outputdir" ]; then mkdir $step06_outputdir; fi

################# Merge the Omni Samples if there are any ##################
omni_samples=false
if [ -e ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes_bplink.bed ]; then
   omni_samples=true
   echo "Merging WGS and Omni genotypes on common markers"
   if [ ! -e merge_bplink_files_v3.sh ]; then ln -s ~/scripts/merge_bplink_files_v3.sh; fi
   bash merge_bplink_files_v3.sh ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes_bplink ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes_bplink ${step06_outputdir}/WGS_Omni_${tissue_of_interest}_merge temp_${tissue_of_interest}_${tmp}
fi

genofile=""
if [ $omni_samples = true ]; then
   genofile="${step06_outputdir}/WGS_Omni_${tissue_of_interest}_merge"
else
   genofile="${step04_outputdir}/${tissue_of_interest}_WGS_genotypes_bplink"
fi

################ Prepare the file for PCA ########################
echo "Running steps before PCA (ie. prePCA.sh)"
if [ ! -e prePCA.sh ]; then ln -s ~/scripts/prePCA.sh; fi
bash prePCA.sh ${genofile} ${step06_outputdir}/${tissue_of_interest}_genofile_prePCA


################# Generate kinship matrix for PCAiR ###############
king2.0 -b ${step06_outputdir}/${tissue_of_interest}_genofile_prePCA_pruned.bed --kinship --prefix ${step06_outputdir}/${tissue_of_interest}_kinship

# Run PC-AiR:
if [ ! -e PC-AiR_v1_empty_kinFile.R ]; then ln -s ~/R_code/PC-AiR_v1_empty_kinFile.R; fi
Rscript PC-AiR_v1_empty_kinFile.R ${step06_outputdir}/${tissue_of_interest}_genofile_prePCA_pruned ${step06_outputdir}/${tissue_of_interest}_kinship ${step06_outputdir}/${tissue_of_interest}_PC-AiR

