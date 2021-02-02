#!/bin/bash -l
#SBATCH -A b2010060
#SBATCH -p core -n 1
#SBATCH -t 05:00:00
#SBATCH -C usage_mail

### Get a table with the gene ID, number of branch-specific non-synonymous substitutions and of synonymous substitutions from PAML results ###
# input files: [gene].fa_rm0_eo_Model2_NSs0_A.results.txt
results="/proj/b2013182/nobackup/bams/Cmon_Ref/PAML/Ns/genes_80_7sp/results" # insert path to directory with results files
mktdir="/proj/b2013182/nobackup/POPseq/mk_tests" # insert path for directory for McDonald-Kreitman tests

cd $mktdir
touch mkt_table_N_S_80_7sp.txt
echo 'gene N*dN S*dS' >> mkt_table_N_S_80_7sp.txt
cd $results
for f in $(ls *.fa_rm0_eo_Model2_NSs0_A.results.txt | sed 's/.fa_rm0_eo_Model2_NSs0_A.results.txt//g'); do
nr=`grep '^TREE #' ${f}.fa_rm0_eo_Model2_NSs0_A.results.txt | awk '{ print $5 }'` # get the ID for C. moneduloides in the gene tree (check for each new PAML run if the position is still correct)
cmo=`echo ${nr:0:1}`
tab=`grep '^  [0-9][0-9]\.\.'$cmo' ' ${f}.fa_rm0_eo_Model2_NSs0_A.results.txt | sed -e 's/^[[:space:]]*//' | awk '{print $8, $9}'` # get N*dN and S*dS from table for branch leading to C. moneduloides
if [[ $tab != "" ]]
then
echo $f $tab >> $mktdir/mkt_table_N_S_80_7sp.txt
fi
done

