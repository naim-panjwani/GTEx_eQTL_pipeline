#!/bin/bash

source 00-set_variables.sh


echo -e "\nStep 02 - Subsetting the RNAseq file for ${tissue_of_interest} samples\n"
if [ ! -d "$step02_outputdir" ]; then mkdir $step02_outputdir; fi
#cat ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt ${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt >${step02_outputdir}/genotyped_RNA_individuals.txt

# Include WGS samples only (if you want Omni samples as well, uncomment out the above and comment out the below):
cat ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt >${step02_outputdir}/genotyped_RNA_individuals.txt

paste <(cat ${step01_outputdir}/${tissue_of_interest}_sampleIDs.txt |cut -d'-' -f1,2) <(cat ${step01_outputdir}/${tissue_of_interest}_sampleIDs.txt |cut -f1) >${step02_outputdir}/sample_individual_IDs_mapping_file.txt
join <(sort ${step02_outputdir}/sample_individual_IDs_mapping_file.txt) <(sort ${step02_outputdir}/genotyped_RNA_individuals.txt) |tr ' ' '\t' |cut -f2 >${step02_outputdir}/genotyped_RNA_samples.txt

(echo -e "Name\nDescription"; cat ${step02_outputdir}/genotyped_RNA_samples.txt) >${step02_outputdir}/gct_columns_to_extract_${tmp}
if [ ! -e subset_gct.py ]; then  ln -s ~/scripts/subset_gct.py; fi
python subset_gct.py ${ExpressionDir}/${ExpressionFile} \
		     ${ExpressionDir}/${ExpressionFileGeneCounts} \
	             ${step02_outputdir}/gct_columns_to_extract_${tmp} \
		     ${step02_outputdir}/${ExpressionFile%.gct.gz}_${tissue_of_interest}


rm ${step02_outputdir}/*${tmp}


