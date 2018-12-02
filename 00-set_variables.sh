tissue_of_interest_search_string="Pancreas"


tissue_of_interest=$(echo ${tissue_of_interest_search_string} |sed 's/ //g' | sed 's/\-/_/g')

tmp="181126"

PhenoDir="../PhenotypeFiles/"
ExpressionDir="../ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/"
WGSDir="../GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1"
GenoArrayDir="../GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/"

SampleAttributesFile="phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz"
#ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v7.p2/pheno_variable_summaries/phs000424.v7.pht002743.v7.GTEx_Sample_Attributes.data_dict.xml
SubjectPhenotypesFile="phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"
#ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v7.p2/pheno_variable_summaries/phs000424.v7.pht002742.v7.GTEx_Subject_Phenotypes.data_dict.xml
SampleAttr_var_ids_file=""
SubjAttr_var_ids_file=""

ExpressionFile="GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
ExpressionFileGeneCounts="GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_reads.gct.gz"
WGSFile="GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC"  # expecting both *.vcf.gz and *_bplink.{bed,bim,fam} files present
GenoArrayFile="GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs" # .vcf.gz
GenoArrayFile_chrX="GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chrX_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs" #.vcf.gz

GENCODE_annotationFile="../gencode.v19.annotation.gtf.gz"

OutputDir="output/"
step01_outputdir="${OutputDir}/01-Intersect_Samples"
step02_outputdir="${OutputDir}/02-RNAseq_${tissue_of_interest}_subset"
step03_outputdir="${OutputDir}/03-Normalized_expression"
step04_outputdir="${OutputDir}/04-subset_WGS_and_Omni"
step05_outputdir="${OutputDir}/05-Run_PEER"
step06_outputdir="${OutputDir}/06-PCA"
step07_outputdir="${OutputDir}/07-Covariates"

includePCA=true
includePEER=true
includeOmni=false
includeSex=true
includeAge=false

