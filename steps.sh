#!/bin/bash

tmp=181130

# Build new 00-set_variables.sh per tissue:
while read -r tissue; do
  echo "tissue_of_interest_search_string=\"${tissue}\"" >tmp_${tmp}
  cp 00-set_variables.sh 00-set_variables.sh.bkp
  cat tmp_${tmp} <(grep -Ev "^tissue_of_interest_search_string" 00-set_variables.sh.bkp) >00-set_variables.sh
  echo ""
  bash 01-check_tissue.sh
  echo ""
done <tissues_of_interest.txt

# After above tissues finish correctly, proceed to step 2 to subset the gene expression matrix for the tissues we want:
# (nohup this in a temporary bash script file
while read -r tissue; do
  echo ""
  echo "tissue_of_interest_search_string=\"${tissue}\"" >tmp_${tmp}
  cp 00-set_variables.sh 00-set_variables.sh.bkp
  cat tmp_${tmp} <(grep -Ev "^tissue_of_interest_search_string" 00-set_variables.sh.bkp) >00-set_variables.sh
  echo ""
  bash 02-subset_gene_expression_for_TOI.sh
  echo ""
done <tissues_of_interest.txt
rm *${tmp}


# Step 3 - Normalizing expression
while read -r tissue; do
  echo ""
  echo "tissue_of_interest_search_string=\"${tissue}\"" >tmp_${tmp}
  cp 00-set_variables.sh 00-set_variables.sh.bkp
  cat tmp_${tmp} <(grep -Ev "^tissue_of_interest_search_string" 00-set_variables.sh.bkp) >00-set_variables.sh
  echo ""
  bash 03-normalize_expression_for_TOI.sh
  echo ""
done <tissues_of_interest.txt

