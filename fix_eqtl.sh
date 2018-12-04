#!/bin/bash

file=$1

paste <(echo -e "chr\tpos"; gunzip -c $file |cut -d' ' -f2 |tr '_' '\t' |awk '{print $1"\t"$2}') <(echo -e "gene_id\tvariant_id\ttss_distance\tpval_nominal\tslope"; gunzip -c $file |tr ' ' '\t') |bgzip -c >${file%.txt.gz}_fixed.txt.gz


