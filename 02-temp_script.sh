tmp=181130

while read -r tissue; do
  echo ""
  echo "tissue_of_interest_search_string=\"${tissue}\"" >tmp_${tmp}
  cp 00-set_variables.sh 00-set_variables.sh.bkp
  cat tmp_${tmp} <(grep -Ev "^tissue_of_interest_search_string" 00-set_variables.sh.bkp) >00-set_variables.sh
  echo ""
  bash 02-subset_gene_expression_for_TOI.sh
  echo ""
done <tissues_of_interest.txt

