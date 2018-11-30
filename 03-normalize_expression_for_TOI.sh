#!/bin/bash

source 00-set_variables.sh
step01_outputdir="${OutputDir}/01-Intersect_Samples"
step02_outputdir="${OutputDir}/02-RNAseq_${tissue_of_interest}_subset"

step03_outputdir="${OutputDir}/03-Normalized_expression"
if [ ! -d "$step03_outputdir" ]; then mkdir $step03_outputdir; fi

if [ ! -e normalize_expression.py ]; then ln -s ~/scripts/normalize_expression.py; fi
python normalize_expression.py ${ExpressionDir}/${ExpressionFile} \
			       ${ExpressionDir}/${ExpressionFileGeneCounts} \
			       ${GENCODE_annotationFile} \
			       ${GenoArrayDir}/${GenoArrayFile}.vcf.gz \
			       ${ExpressionFile%.gct.gz}_normalized \
			       --output_dir=$step03_outputdir \
			       --expression_threshold=0.1 \
			       --count_threshold=5 \
			       --min_samples=10
if [ ! -e fix_bed.py]; then ln -s ~/scripts/fix_bed.py; fi
python fix_bed.py ${step03_outputdir}/${ExpressionFile%.gct.gz}_normalized.bed.gz

# Replace fixed bed and remove the old:
rm ${step03_outputdir}/${ExpressionFile%.gct.gz}_normalized.bed.gz
mv ${ExpressionFile%.gct.gz}_normalized_fixed.bed.gz ${ExpressionFile%.gct.gz}_normalized.bed.gz

