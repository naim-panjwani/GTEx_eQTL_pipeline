#!/bin/bash

# Symbolic link to GTEx V8 dbGaP folder:
ln -s /hpf/largeprojects/struglis/datasets/Public/GTEx/GTEx_V8_dbGaP/

# Create local data folders and symbolic links
mkdir data
mkdir data/genotypes
mkdir data/expression
mkdir data/phenotypes
mkdir data/phenotypes/sample_subsets
mkdir data/expression/bedfile_subsets
mkdir output/
ln -s /hpf/largeprojects/struglis/datasets/Public/GTEx/GTEx_V8_dbGaP/Phenotype_Files/phs000424.v8.pht002742.v8.p2.GTEx_Subject_Phenotypes.var_report.xml data/phenotypes/GTEx_Subject_Phenotypes.var_report.xml
ln -s /hpf/largeprojects/struglis/datasets/Public/GTEx/GTEx_V8_dbGaP/Phenotype_Files/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz data/phenotypes/GTEx_Subject_Phenotypes.GRU.txt.gz
ln -s /hpf/largeprojects/struglis/datasets/Public/GTEx/GTEx_V8_dbGaP/Phenotype_Files/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz data/phenotypes/GTEx_Sample_Attributes.GRU.txt.gz
ln -s /hpf/largeprojects/struglis/datasets/Public/GTEx/GTEx_V8_dbGaP/Phenotype_Files/phs000424.v8.pht002743.v8.GTEx_Sample_Attributes.data_dict.xml data/phenotypes/GTEx_Sample_Attributes.data_dict.xml

# Extract the necessary files that will be needed:
tar -xvf GTEx_V8_dbGaP/Expression_Files/phe000037.v1.GTEx_v8_RNAseq.expression-data-matrixfmt.c1.GRU.tar -C data/expression/
tar -xvf GTEx_V8_dbGaP/Genotype_Files/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar -C data/genotypes/


# Download files needed:
## Gene read counts
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz -P data/expression/
## Gene TPM values
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz -P data/expression/
## Gene annotations (collapsed):
wget https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf -P data/expression/

# Convert vcf file to PLINK
vcf="data/genotypes/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
plink2 --vcf "$vcf" --make-bed --out "${vcf%.vcf.gz}_bplink" --maf 0.01 --memory 10000 --threads 2


# Make covariates file
## Get sequencing platform (HiSeq 2000 vs HiSeq X) and sequencing protocol (PCR-based vs PCR-free) per Subject
gunzip -c data/phenotypes//GTEx_Sample_Attributes.GRU.txt.gz |grep -v ^"#" |cut -f2,25,26 |grep -E "DNA" |grep -E "WGS" |cut -f1,2 > tmp
paste <(cut -f1 tmp |tr '-' '\t' |cut -f1-2 |tr '\t' '-') <(cut -f2 tmp |cut -d' ' -f1) <(cut -f2 tmp |cut -d' ' -f2- |sed 's/v1\|v2//g' |tr -s " ") > data/phenotypes/GTEx_Subjects_sequencing_covariates.txt
(echo -e "SUBJID\tProtocol\tPlatform"; cat data/phenotypes/GTEx_Subjects_sequencing_covariates.txt) > tmp
cp tmp data/phenotypes/GTEx_Subjects_sequencing_covariates.txt
rm tmp

## Add Sex and Age covariates
(echo -e "SUBJID\tProtocol\tPlatform\tSEX\tAGE"; join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,2.2,2.3 <(sed '1d' data/phenotypes/GTEx_Subjects_sequencing_covariates.txt |sort -k1)  <(gunzip -c data/phenotypes/GTEx_Subject_Phenotypes.GRU.txt.gz |grep -v ^"#" |sed '1,2d' |cut -f2,4-5 |sort -k1)) > data/phenotypes/GTEx_Subjects_covariates.txt
sed 's/\ /_/g' data/phenotypes/GTEx_Subjects_covariates.txt |sed 's/(//g' |sed 's/)//g' > tmp
cp tmp data/phenotypes/GTEx_Subjects_covariates.txt
rm tmp

## transpose file
bash ~/scripts/transposing/transpose_with_awk.sh data/phenotypes/GTEx_Subjects_covariates.txt tmp
sed -i 's/SUBJID/ID/g' tmp
cp tmp data/phenotypes/GTEx_Subjects_covariates.txt
rm tmp

# Use R to convert the factors to integers:
