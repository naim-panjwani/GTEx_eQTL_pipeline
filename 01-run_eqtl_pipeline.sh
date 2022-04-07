#!/bin/bash

export PATH="/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv2/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin:${PATH}"
#source activate rbioenv2

sample_participant_lookup=$1
num_threads=$2

i="$sample_participant_lookup"
rand=$RANDOM

vcf="data/genotypes/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
vcf_chr_list="data/genotypes/chromlist.txt"
counts_gct="data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
tpm_gct="data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
annotation_gtf="data/expression/gencode.v26.GRCh38.genes.gtf"


# Generate tabix index file for $vcf
[ ! -f "${vcf}.tbi" ] && tabix -p vcf "$vcf"
# Generate list of chromosomes in $vcf
[ ! -f "${vcf_chr_list}" ] && tabix -l "$vcf" > "$vcf_chr_list"


# Get GTEx V8 pipeline docker image:
module load Singularity
[ ! -f gtex_eqtl_V8.sif ] && singularity pull docker://broadinstitute/gtex_eqtl:V8 # a fix was required on read_gct() to set the index column to 'Name' instead of 0. Use gtex_eqtl_V8_fixed.sif instead

# For each sample subset:
#sample_participant_lookup="data/phenotypes/sample_subsets/Brain.Cortex.female.SAMPID_lookup.txt"
#for i in $(ls data/phenotypes/sample_subsets/*SAMPID_lookup.txt); do
    echo ""
    echo "Processing eQTL analysis for $i"
    echo ""
    #sample_participant_lookup="$i"
    tissue=$(basename ${sample_participant_lookup} |sed -E 's/\.SAMPID_list.txt|\.SAMPID_lookup.txt|\.sample_attributes.txt//g' |sed -E 's/\.female|\.male//g' |sed 's/\./_/g' |tr A-Z a-z)
    tissue=$(echo "$tissue" |sed 's/brain_spinal_cord_cervical_c1/brain_spinal_cord_cervical_c-1/g')
    # Download gene read counts and TPMs for $tissue
    counts_gct="data/expression/bedfile_subsets/gene_reads_2017-06-05_v8_${tissue}.gct.gz"
    tpm_gct="data/expression/bedfile_subsets/gene_tpm_2017-06-05_v8_${tissue}.gct.gz"
#     if [ ! -f "${counts_gct}" ]; then
#         wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_${tissue}.gct.gz -P data/expression/bedfile_subsets/
#         wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_${tissue}.gct.gz -P data/expression/bedfile_subsets/
#     fi

    # Generate BED expression files
    prefix=$(basename $sample_participant_lookup |sed 's/.SAMPID_lookup.txt//g')
    [ ! -d "output/${prefix}" ] && mkdir "output/${prefix}"
    prefix="output/${prefix}/${prefix}"
    singularity run gtex_eqtl_V8_fixed.sif /bin/bash -c "/src/eqtl_prepare_expression.py \
        ${tpm_gct} ${counts_gct} ${annotation_gtf} \
        ${sample_participant_lookup} ${vcf_chr_list} ${prefix} \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm"

    # Calculate PEER factors
    sample_size=$(gunzip -c ${prefix}.expression.bed.gz |head -1 |wc -w)
    sample_size=$(( sample_size - 4 ))
    num_peer=0
    if [[ $sample_size -lt 150 ]]; then
        num_peer=15
    elif [[ $sample_size -ge 150 ]] && [[ $sample_size -lt 250 ]]; then
        num_peer=30
    elif [[ $sample_size -ge 250 ]] && [[ $sample_size -lt 350 ]]; then
        num_peer=45
    elif [[ $sample_size -ge 350 ]]; then
        num_peer=60
    else
        echo "Invalid sample size: ${sample_size}"
    fi
    echo "Sample size: ${sample_size}"
    echo "Number of PEER factors: ${num_peer}"
    singularity run gtex_eqtl_V8_fixed.sif /bin/bash -c "Rscript /src/run_PEER.R ${prefix}.expression.bed.gz ${prefix} ${num_peer}"


    # Calculate genotype PC's
    vcf_bplink="${vcf%.vcf.gz}_bplink"
    sample_file="${prefix}_samples.txt"
    head -1 "${prefix}.PEER_covariates.txt" |tr '\t' '\n' |sed '1d' |awk '{print 0"\t"$1}' > "$sample_file"
    bash ~/scripts/doPCA.sh "$sample_file" "$vcf_bplink" "${prefix}_PCAiR"
    echo -e "$(basename $prefix)\t$(wc -l ${prefix}_PCAiR_eigenvectors.txt |cut -d' ' -f1)" >> output/tissue_sample_sizes.txt
    # invert eigenvectors file (top 5 PC's only):
    bash ~/scripts/transposing/transpose_with_awk.sh <(cut -d' ' -f1-6 ${prefix}_PCAiR_eigenvectors.txt) "tmp${rand}"
    paste <(echo -e "ID\nPC1\nPC2\nPC3\nPC4\nPC5") "tmp${rand}" > "${prefix}.PCAiR_covariates.txt"
    rm "tmp${rand}"

    # Combine all covariates
    covarfile1="${prefix}.PEER_covariates.txt"
    covarfile2="${prefix}.PCAiR_covariates.txt"
    covarfile3="data/phenotypes/GTEx_Subjects_covariates_numeric.txt"

    singularity run gtex_eqtl_V8_fixed.sif /bin/bash -c "/src/combine_covariates.py \
        $covarfile1 $prefix \
        --genotype_pcs ${covarfile2} \
        --add_covariates ${covarfile3}"

    singularity run gtex_eqtl_V8_fixed.sif /bin/bash -c "/opt/fastqtl/python/run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} \
        --covariates ${prefix}.combined_covariates.txt \
        --window 1e6 --chunks 100 --threads ${num_threads} \
        --maf_threshold 0.01 --ma_sample_threshold 5"
#done
