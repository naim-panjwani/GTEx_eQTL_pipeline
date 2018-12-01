#!/bin/bash

tmp=181201

source 00-set_variables.sh
if [ ! -d "$step04_outputdir" ]; then mkdir $step04_outputdir; fi

# Make the keep file for plink:
paste ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt >${step04_outputdir}/bplink_keep_file_${tissue_of_interest}_${tmp}
paste ${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt ${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt >${step04_outputdir}/bplink_keep_file_Omni_${tissue_of_interest}_${tmp}

# Extract with plink:
plink --bfile ${WGSDir}/${WGSFile}_bplink --keep ${step04_outputdir}/bplink_keep_file_${tissue_of_interest}_${tmp} --make-bed --out ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes_bplink
plink --bfile ${GenoArrayDir}/${GenoArrayFile}_bplink_iid_update --keep ${step04_outputdir}/bplink_keep_file_Omni_${tissue_of_interest}_${tmp} --make-bed --out ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes_bplink

# Convert to vcf:
plink --bfile ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes_bplink --recode vcf --out ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes
plink --bfile ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes_bplink --recode vcf --out ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes

# Compress & index:
bgzip ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes.vcf
tabix -p vcf ${step04_outputdir}/${tissue_of_interest}_WGS_genotypes.vcf.gz
bgzip ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes.vcf
tabix -p vcf ${step04_outputdir}/${tissue_of_interest}_Omni_genotypes.vcf.gz


rm ${step04_outputdir}/*${tmp}

