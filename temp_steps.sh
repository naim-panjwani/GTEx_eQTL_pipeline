#!/bin/bash

# Genotypes:
ln -s ../../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
ln -s ../../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz.tbi

# Phenotypes (rank normalized RPKM)
ln -s ../../ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/normalized_rpkm_output/GTEx_v7_rpkm_normalized_Pancreas.expression.bed.gz

# Covariates:
mkdir covariates_building
ln -s ../../../ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/normalized_rpkm_output/GTEx_v7_rpkm_normalized_Pancreas.PEER_covariates.txt
ln -s ../../../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_bplink_PC-AiR_eigenvectors.txt
ln -s ../../PhenotypeFiles/Samples_WGS_attributes.txt

# Creating the covariates file in R:
mkdir covariates_building
cd covariates_building/
./create_covars.R
bgzip Covariates_with_age.txt
bgzip Covariates_without_age.txt

# Subset and order the genotype and phenotype (expression) files to match:
gunzip -c GTEx_v7_rpkm_normalized_Pancreas.expression.bed.gz |head -1 |cut -f5- |tr '\t' '\n' >Panc_RNAseq_subjects.txt # header line = 1
headerLine=$(gunzip -c GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz |head -300 |grep -n ^"#CHROM" |cut -f1 |cut -d':' -f1) # 217
gunzip -c GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz |head -${headerLine} |tail -1 |cut -f10- |tr '\t' '\n' >WGS_subjects.txt

#(gunzip -c GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz |head -${headerLine} |tail -1 |cut -f1-9 |tr '\t' '\n'; comm -12 <(sort WGS_subjects.txt) <(sort Panc_RNAseq_subjects.txt)) >Genotype_columns_to_cut.txt

(gunzip -c GTEx_v7_rpkm_normalized_Pancreas.expression.bed.gz |head -1 |cut -f1-4 |tr '\t' '\n'; comm -12 <(sort Panc_RNAseq_subjects.txt) <(sort WGS_subjects.txt)) >Pancreas_expression_columns_to_cut.txt

###########################################################################################
#######    The VCF file does not need to be broken down:    ###############################
###########################################################################################
# best to subset the region first as the file is huge (120GB)
#tabix -h GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz X:110000000-120000000 |bgzip -c >GTEx_WGS_PLINKQC_X_110-120Mb.vcf.gz
#python ~/scripts/cut2.py -i GTEx_WGS_PLINKQC_X_110-120Mb.vcf.gz \
#			 -o GTEx_WGS_PLINKQC_X_110-120Mb_Pancreas_subjects.vcf.gz \
#			--columnsfile=Genotype_columns_to_cut.txt \
#			--headerLine=${headerLine}
# Takes tooo long to try to do the entire VCF file (I killed it after two days of running)

python ~/scripts/cut2.py -i GTEx_v7_rpkm_normalized_Pancreas.expression.bed.gz \
			 -o Pancreas_expression.bed.gz \
			--columnsfile=Pancreas_expression_columns_to_cut.txt \
			--headerLine=1


#gunzip -c Pancreas_expression.bed.gz |bgzip -c >Pancreas_expression2.bed.gz
#rm Pancreas_expression.bed.gz
#mv Pancreas_expression2.bed.gz Pancreas_expression.bed.gz
#tabix -p bed Pancreas_expression.bed.gz

# Waiting for cut2.py to finish
gunzip -c Pancreas_genotypes.vcf.gz |bgzip -c >Pancreas_genotypes2.vcf.gz
rm Pancreas_genotypes.vcf.gz
mv Pancreas_genotypes2.vcf.gz Pancreas_genotypes.vcf.gz
tabix -p vcf Pancreas_genotypes.vcf.gz

# Filter out genotypes with MAF<1%
bplink="../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_bplink"
cat Genotype_columns_to_cut.txt |sed '1,9d' |paste - <(cat Genotype_columns_to_cut.txt |sed '1,9d') >Genotype_columns_to_cut2.txt
plink1.9 --bfile $bplink --chr 23 --from-kb 114000 --to-kb 117000 --maf 0.01 --make-bed --out chr23_114-117Mb_WGS_Pancreas_samples_maf1pct --keep Genotype_columns_to_cut2.txt
plink1.9 --bfile chr23_114-117Mb_WGS_Pancreas_samples_maf1pct --freq --out chr23_114-117Mb_WGS_Pancreas_samples_freq --keep Genotype_columns_to_cut2.txt
plink1.9 --bfile chr23_114-117Mb_WGS_Pancreas_samples_maf1pct --recode vcf-iid --out chr23_114-117Mb_WGS_Pancreas_samples_maf1pct_vcf
sed 's/^23/X/g' chr23_114-117Mb_WGS_Pancreas_samples_maf1pct_vcf.vcf |bgzip -c >chrX_114-117Mb_WGS_Pancreas_samples_maf1pct.vcf.gz
tabix -p vcf chrX_114-117Mb_WGS_Pancreas_samples_maf1pct.vcf.gz
rm chr23_114-117Mb_WGS_Pancreas_samples_maf1pct_vcf.vcf

# Run FastQTL:
vcffile="chrX_114-117Mb_WGS_Pancreas_samples_maf1pct.vcf.gz"
expressionfile="Pancreas_expression.bed.gz"
chr="X"
startbp=114000000
endbp=117000000
outfile="eQTL_results/${vcffile%.vcf.gz}"
mkdir eQTL_results


cd ../
ln -s ../../FastQTL/bin/fastQTL.static
./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --permute 1000 10000 \
                 --out ${outfile}_permutations.txt.gz \
                 --cov covariates_building/Covariates_without_age.txt.gz \
		 --window 2e6

# The permutation test only outputs the top variant (for the case of SLC6A14, rs378766)
./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --permute 1000 10000 \
                 --out ${outfile}_permutations_age_adjusted.txt.gz \
                 --cov covariates_building/Covariates_with_age.txt.gz \
                 --window 2e6

./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --out ${outfile}_FastQTL.txt.gz \
                 --cov covariates_building/Covariates_without_age.txt.gz \
                 --window 2e6

./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --out ${outfile}_FastQTL_age_adjusted.txt.gz \
                 --cov covariates_building/Covariates_with_age.txt.gz \
                 --window 2e6
# Add chr pos columns, make it tab-separated, and add standard error column:
cd eQTL_results
for i in $(ls *txt.gz |grep -Ev "*fixed*|*permutations*"); do bash fix_bed.sh $i; done
for i in $(ls *txt.gz |grep -E "*fixed*"); do Rscript include_stderr.R $i; done
for i in $(ls *slopeSE.txt); do bgzip $i; done
for i in $(ls *txt.gz |grep -Ev "*fixed_slopeSE*|*permutations*"); do rm $i; done


# Incorporate the minor allele count and frequency for the RNAseq samples analyzed per variant:
# Will do this later...


# Running SLC6A14 Pancreas by sex:
./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --out ${outfile}_FastQTL_males.txt.gz \
                 --cov covariates_building/Covariates_without_age.txt.gz \
                 --window 2e6 \
		 --include-samples Pancreas_male_samples2.txt


./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --out ${outfile}_FastQTL_females.txt.gz \
                 --cov covariates_building/Covariates_without_age.txt.gz \
                 --window 2e6 \
                 --include-samples Pancreas_female_samples2.txt

cd eQTL_results/
bash fix_bed.sh chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_males.txt.gz
bash fix_bed.sh chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females.txt.gz
Rscript include_stderr.R chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_males_fixed.txt.gz
Rscript include_stderr.R chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_fixed.txt.gz
bgzip chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_males_fixed_slopeSE.txt
bgzip chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_fixed_slopeSE.txt
rm chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_males_fixed.txt.gz
rm chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_males.txt.gz 
rm chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_fixed.txt.gz 
rm chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females.txt.gz


# What is the genotype distribution of these individuals at rs3788766 - males and females separated
plink1.9 --bfile ../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_bplink --keep Pancreas_male_samples3.txt --snp X_115566839_A_G_b37 --recodeA --out Pancreas_male_samples_rs3788766_genotypes
plink1.9 --bfile ../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_bplink --keep Pancreas_female_samples3.txt --snp X_115566839_A_G_b37 --recodeA --out Pancreas_female_samples_rs3788766_genotypes

cut -d' ' -f7 Pancreas_male_samples_rs3788766_genotypes.raw |sed '1d' |sort |uniq -c
#     75 0
#     54 2
#      3 NA

cut -d' ' -f7 Pancreas_female_samples_rs3788766_genotypes.raw |sed '1d' |sort |uniq -c
#     31 0
#     45 1
#     18 2


# Conduct eQTL analysis on SLC6A14 region on females without those who are heterozygous at rs3788766:
./fastQTL.static --vcf $vcffile \
                 --bed $expressionfile \
                 --region ${chr}:${startbp}-${endbp} \
                 --out ${outfile}_FastQTL_females_rs3788766_hets_rm.txt.gz \
                 --cov covariates_building/Covariates_without_age.txt.gz \
                 --window 2e6 \
                 --include-samples Pancreas_female_samples_rs3788766_hets_rm.txt
cd eQTL_results/
bash fix_bed.sh chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_rs3788766_hets_rm.txt.gz
Rscript include_stderr.R chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_rs3788766_hets_rm_fixed.txt.gz
bgzip chrX_114-117Mb_WGS_Pancreas_samples_maf1pct_FastQTL_females_rs3788766_hets_rm_fixed_slopeSE.txt








