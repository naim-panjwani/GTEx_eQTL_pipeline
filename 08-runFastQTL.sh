#!/bin/bash

tmp=181201
source 00-set_variables.sh
if [ ! -d "$step08_outputdir" ]; then mkdir $step08_outputdir; fi

chr=$1
startbp=$2
endbp=$3
covariatesfile=$4
outsuffix=$5
if [ -z $covariatesfile ]; then covariatesfile="${step07_outputdir}/${tissue_of_interest}_covariates.txt"; fi
if [ ! -z "$outsuffix" ]; then outsuffix="_${outsuffix}"; fi
vcffile="${step04_outputdir}/${tissue_of_interest}_WGS_genotypes"

expressionfile="${step03_outputdir}/${ExpressionFile%.gct.gz}_${tissue_of_interest}_normalized.expression.bed.gz"

# Filter out variants with MAF<1% and isolate the region of interest:
pchr=$chr
if [ $chr = "X" ]; then pchr=23; fi
if [ $chr = 23 ]; then chr="X"; pchr=23; fi
plink --bfile ${vcffile}_bplink --maf 0.01 --chr ${pchr} --from-bp ${startbp} --to-bp ${endbp} --make-bed --out ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_${tmp}
#plink --bfile ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp} --freq --out ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_freq

#if [ $chr = 23 ]; then chr="X"; fi
plink --bfile ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_${tmp} --recode vcf-iid --out ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf
if [ $chr = "X" ]; then
  sed 's/^23/X/g' ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf |bgzip -c >${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf_chrFixed.vcf.gz
  rm ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf
  mv ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf_chrFixed.vcf.gz ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf.gz
else
  bgzip ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf
fi
tabix -p vcf ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf.gz

# Run FastQTL:
if [ ! -e fastQTL.static ]; then ln -s ../FastQTL/bin/fastQTL.static; fi
./fastQTL.static --vcf ${step08_outputdir}/${tissue_of_interest}_bplink_maf1pct_${chr}_${startbp}_${endbp}_vcf.vcf.gz \
		 --bed $expressionfile \
		 --region ${chr}:${startbp}-${endbp} \
		 --cov $covariatesfile \
		 --window 2e6 \
		 --out ${step08_outputdir}/${tissue_of_interest}_FastQTL_${chr}_${startbp}_${endbp}${outsuffix}.txt
bgzip ${step08_outputdir}/${tissue_of_interest}_FastQTL_${chr}_${startbp}_${endbp}${outsuffix}.txt
bash fix_eqtl.sh ${step08_outputdir}/${tissue_of_interest}_FastQTL_${chr}_${startbp}_${endbp}${outsuffix}.txt.gz
Rscript include_stderr.R ${step08_outputdir}/${tissue_of_interest}_FastQTL_${chr}_${startbp}_${endbp}${outsuffix}_fixed.txt.gz
bgzip ${step08_outputdir}/${tissue_of_interest}_FastQTL_${chr}_${startbp}_${endbp}${outsuffix}_fixed_slopeSE.txt


rm ${step08_outputdir}/${tissue_of_interest}*${tmp}*


