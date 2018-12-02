#!/bin/bash

source 00-set_variables.sh
if [ ! -d "$step05_outputdir" ]; then mkdir $step05_outputdir; fi


# Run PEER factors on the bed file output by normalize_expression.py script:
bedfile="${step03_outputdir}/${ExpressionFile%.gct.gz}_${tissue_of_interest}_normalized.expression.bed.gz"
prefix="${ExpressionFile%.gct.gz}_${tissue_of_interest}_normalized"
num_samples=$(gunzip -c ${bedfile} |head -1 |cut -f5- |tr '\t' '\n' |wc -l)
num_peer=15
if [ $num_samples -ge 350 ]; then 
   num_peer=60
elif [ $num_samples -ge 250 ]; then
   num_peer=45
elif [ $num_samples -ge 150 ]; then
   num_peer=30
else
   num_peer=15
fi

echo "Running PEER for ${prefix}.expression.bed.gz"
if [ ! -e run_PEER.R ]; then ln -s ~/scripts/GTEx_scripts/run_PEER.R; fi
Rscript run_PEER.R ${bedfile} ${prefix} ${num_peer} -o ${step05_outputdir}

