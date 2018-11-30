#!/bin/bash

source 00-set_variables.sh

# File/directory checks (work in progress)
if [ ! -d "${OutputDir}/" ]; then mkdir ${OutputDir}/; fi

# Gather sample ID's for tissue of interest
echo "Step 01 - Intersecting RNAseq, WGS, and Omni Sample/Individual IDs"
step01_outputdir="${OutputDir}/01-Intersect_Samples"
if [ ! -d "$step01_outputdir" ]; then mkdir $step01_outputdir; fi
#echo "Extracting Sample IDs from sample attributes file for ${tissue_of_interest} to ${step01_outputdir}/${tissue_of_interest}_sampleIDs.txt"
gunzip -c "${PhenoDir}/${PhenoFile}"  |cut -f2,13-15,17 >${step01_outputdir}/SampleID_and_tissue_Attributes_${tmp}
grep "${tissue_of_interest_search_string}" ${step01_outputdir}/SampleID_and_tissue_Attributes_${tmp} >${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt
num_tissue_types=$(cut -f4 ${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt |sort |uniq -c |wc -l)
if [ $num_tissue_types -gt 1 ]; then 
   echo "Have more than one tissue type!! Check the ${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt file"
   echo "Change search string to one of:"
   cut -f4 ${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt |sort |uniq
   exit 1
fi

num_tissue_samples=$(wc -l ${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt |cut -d' ' -f1)
gunzip -c ${ExpressionDir}/${ExpressionFile} |sed '1,2d' |head -1 |cut -f3- |tr '\t' '\n' >${step01_outputdir}/RNAseq_sampleIDs.txt
comm -12 <(cut -f1 ${step01_outputdir}/${tissue_of_interest}_sampleID_attributes.txt |sort) <(sort ${step01_outputdir}/RNAseq_sampleIDs.txt) >${step01_outputdir}/${tissue_of_interest}_sampleIDs.txt
echo "Number of ${tissue_of_interest} samples collected: ${num_tissue_samples}"


#echo "Extracting Individual ID's from Sample ID's to ${step01_outputdir}/${tissue_of_interest}_individual_IDs.txt"
cat ${step01_outputdir}/${tissue_of_interest}_sampleIDs.txt |cut -d'-' -f1,2 >${step01_outputdir}/${tissue_of_interest}_individual_IDs.txt
num_with_RNAseq="$(wc -l ${step01_outputdir}/${tissue_of_interest}_individual_IDs.txt |cut -d' ' -f1)"
echo "Total number of individuals with RNAseq: ${num_with_RNAseq}"


#echo "Listing genotyped individuals in WGS_IDs.txt and imputed_samples.txt"
gunzip -c ${WGSDir}/${WGSFile}.vcf.gz |head -300 |grep ^"#CHROM" |cut -f10- |tr '\t' '\n' >${step01_outputdir}/WGS_IDs_${tmp}
gunzip -c ${GenoArrayDir}/${GenoArrayFile}.vcf.gz |head -300 |grep ^"#CHROM" |cut -f10- |tr '\t' '\n' >${step01_outputdir}/imputed_samples_${tmp}
cut -d'-' -f1,2 ${step01_outputdir}/imputed_samples_${tmp} >${step01_outputdir}/imputed_samples2_${tmp}


#echo "Listing the 34 individuals with only Omni genotype data"
#echo "Note: There are actually 42 Omni genotyped individuals,"
#echo "but 8 are excluded as identified in the WGS for having large chromosomal abnormalities"
#echo ""
#echo "List of 34 individuals with only Omni data is provided in imputed_samples_not_WGS.txt"
grep "EXCLUDE" ${PhenoDir}/Samples_WGS_attributes.txt |cut -f1 |cut -d'-' -f1,2 >${step01_outputdir}/WGS_excluded_IDs_${tmp}  # ID's excluded in WGS due to large CNVs, trisomies, and mosaicisms
# We do not want the genotypes of the 17 individuals with large genomic structural problems, so let's exclude these genotypes:
comm <(sort ${step01_outputdir}/imputed_samples2_${tmp}) <(sort ${step01_outputdir}/WGS_excluded_IDs_${tmp}) |cut -f1 |awk 'NF>0' >${step01_outputdir}/imputed_samples3_${tmp}
# List only the ones that do not have WGS done:
comm <(sort ${step01_outputdir}/imputed_samples3_${tmp}) <(sort ${step01_outputdir}/WGS_IDs_${tmp}) |cut -f1 |awk 'NF>0' >${step01_outputdir}/imputed_samples_not_WGS_${tmp} # and 33/34 of these have RNAseq btw


#echo "Listing samples with both WGS and RNAseq in ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt"
comm -12 <(cat ${step01_outputdir}/${tissue_of_interest}_individual_IDs.txt |cut -d'-' -f1,2 |sort) <(sort  ${step01_outputdir}/WGS_IDs_${tmp}) >${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt
num_with_WGS=$(wc -l ${step01_outputdir}/${tissue_of_interest}_individual_IDs_intersection_WGS.txt |cut -d' ' -f1)
echo "${num_with_WGS} RNAseq samples found in WGS dataset"


#echo "Listing samples for which we only have genotypes from microarray in ${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt"
comm -12 <(cat ${step01_outputdir}/${tissue_of_interest}_individual_IDs.txt |cut -d'-' -f1,2 |sort) <(sort ${step01_outputdir}/imputed_samples_not_WGS_${tmp}) >${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt
num_Omni=$(wc -l ${step01_outputdir}/${tissue_of_interest}_individual_IDs_in_Omni_only.txt |cut -d' ' -f1)
echo "${num_Omni} RNAseq samples found in Omni dataset"

echo "$(( $num_with_WGS + $num_Omni )) available for eQTL analysis"

rm ${step01_outputdir}/*${tmp}

