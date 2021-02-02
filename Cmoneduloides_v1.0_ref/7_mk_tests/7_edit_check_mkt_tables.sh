#!/bin/bash

mktdir="/proj/b2013182/nobackup/POPseq/mk_tests"
cd $mktdir

# remove " from tables
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*); do sed -i.bak 's/"//g' ${i}; done

# compare missing gene lists from R output to original divergence tables before merging with the polymorphism table
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_N_S.txt`; do grep "^$i " mkt_table_N_S_80_5sp_sorted_unique.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_N_S_doublecheck.txt; done
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_N_S.txt`; do grep "^$i " mkt_table_N_S_80_7sp_sorted_unique.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_N_S_doublecheck.txt; done
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_N_S.txt`; do grep "^$i " mkt_table_N_S_80_8sp_H_NC_sorted_unique.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_N_S_doublecheck.txt; done

# compare missing gene lists from R output to results files in PAML directories
cd /proj/b2013182/nobackup/bams/Cmon_Ref/PAML/Ns/genes_80_5sp/results
for i in `awk '{print $1}' $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_N_S.txt`; do ls ${i}.fa_rm0_eo_Model2_NSs0_A.results.txt >> $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_N_S_vs_PAML_Model2_results_files.txt; done

cd /proj/b2013182/nobackup/bams/Cmon_Ref/PAML/Ns/genes_80_7sp/results
for i in `awk '{print $1}' $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_N_S.txt`; do ls ${i}.fa_rm0_eo_Model2_NSs0_A.results.txt >> $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_N_S_vs_PAML_Model2_results_files.txt; done

cd /proj/b2013182/nobackup/bams/Cmon_Ref/PAML/Ns/genes_80_8sp_H_NC/results
for i in `awk '{print $1}' $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_N_S.txt`; do ls ${i}.fa_rm0_eo_Model2_NSs0_A.results.txt >> $mktdir/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_N_S_vs_PAML_Model2_results_files.txt; done

# compare missing gene lists from R output to original polymorphism table before merging with the divergence tables
cd $mktdir
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_pN_pS.txt`; do grep "Cmoneduloides_001-${i}" Cmon_allSites_genes_annotations_160818.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_5sp_sorted_unique_missing_genes_pN_pS_vs_allSites_vcf_annotation_file.txt; done
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_pN_pS.txt`; do grep "Cmoneduloides_001-${i}" Cmon_allSites_genes_annotations_160818.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_pN_pS_vs_allSites_vcf_annotation_file.txt; done
for i in `awk '{print $1}' mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_pN_pS.txt`; do grep "Cmoneduloides_001-${i}" Cmon_allSites_genes_annotations_160818.txt >> mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_missing_genes_pN_pS_vs_allSites_vcf_annotation_file.txt; done

### MANUALLY CHECK EACH OUTPUT FILE FROM THIS SCRIPT BEFORE MOVING ON ###

