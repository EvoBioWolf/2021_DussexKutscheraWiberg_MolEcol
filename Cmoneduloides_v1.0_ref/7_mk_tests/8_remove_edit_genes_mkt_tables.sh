#!/bin/bash

mktdir="/proj/b2013182/nobackup/POPseq/mk_tests"
cd $mktdir

# remove all genes with missing data from the large mkt tables
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*_sorted_unique_180118.txt | sed 's/.txt//g'); do grep -v " NA " ${i}.txt > ${i}_noMissingData.txt; done

# set genes that are present in the VCF file but without variable sites to 0 (excluding multiple stopp codons)
# make a list of genes to be set to 0
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*sorted_unique_missing_genes_pN_pS_vs_allSites_vcf_annotation_file.txt | sed 's/_vs_allSites_vcf_annotation_file.txt//g'); do awk -F'\t' '$18 != "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS" {print $16}' ${i}_vs_allSites_vcf_annotation_file.txt | sed 's/Cmoneduloides_001-//g' | sort -u > ${i}_setTo0_gene_list.txt; done

# set the genes from the lists to 0
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*sorted_unique_missing_genes_pN_pS.txt | sed 's/.txt//g'); do for j in `cat ${i}_setTo0_gene_list.txt`; do grep ${j} ${i}.txt | sed 's/ NA/ 0/g' >> ${i}_setTo0.txt; done; done

# merge the recoded tables with the large mkt tables (noMissingData)
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*_sorted_unique_180118.txt | sed 's/.txt//g'); do cat ${i}_noMissingData.txt ${i}_missing_genes_pN_pS_setTo0.txt | sort -u > ${i}_missingSetTo0.txt; done

# 2018-01-18: no gene found that has to be set to 0, simply renamed the large mkt tables
for i in $(ls mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_*_sorted_unique_180118.txt | sed 's/.txt//g'); do cp ${i}_noMissingData.txt ${i}_missingSetTo0.txt; done