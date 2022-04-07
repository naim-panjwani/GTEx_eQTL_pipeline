#!/bin/bash

[ ! -f output/tissue_sample_sizes.txt ] && touch output/tissue_sample_sizes.txt

for i in $(ls data/phenotypes/sample_subsets/*SAMPID_lookup.txt); do
    echo ""
    echo "Processing eQTL analysis for $i"
    echo ""
    sample_participant_lookup="$i"
    jobname=$(basename $sample_participant_lookup |sed 's/.SAMPID_lookup.txt//g')
    tissue=$(basename ${sample_participant_lookup} |sed -E 's/\.SAMPID_list.txt|\.SAMPID_lookup.txt|\.sample_attributes.txt//g' |sed -E 's/\.female|\.male//g' |sed 's/\./_/g' |tr A-Z a-z)
    tissue=$(echo "$tissue" |sed 's/brain_spinal_cord_cervical_c1/brain_spinal_cord_cervical_c-1/g')

    # Downloads must happen in interactive node before submitting job:
    counts_gct="data/expression/bedfile_subsets/gene_reads_2017-06-05_v8_${tissue}.gct.gz"
    tpm_gct="data/expression/bedfile_subsets/gene_tpm_2017-06-05_v8_${tissue}.gct.gz"
    if [ ! -f "${counts_gct}" ]; then
        wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_${tissue}.gct.gz -P data/expression/bedfile_subsets/
        wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_${tissue}.gct.gz -P data/expression/bedfile_subsets/
    fi
    
    qsub -F "${i} 2" 01-run_eqtl_pipeline.sh -l walltime=23:59:00 -l nodes=1:ppn=4 -l mem=20g -l vmem=20g -o ./jobout -e ./jobout -d `pwd` -N "${jobname}"
done
