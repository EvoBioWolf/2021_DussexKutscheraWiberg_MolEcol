#!/bin/bash

# !!! Before running this script the input table needs to be manually edited according to README file in /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mk_tests/ !!!

cd /proj/b2013182/nobackup/POPseq/mk_tests

#extract gene annotations (maker) from large table, check number of genes
head -n1 Cmon_allSites_annotations_160818.txt > Cmon_allSites_genes_annotations_160818.txt
grep 'Cmoneduloides_001' Cmon_allSites_annotations_160818.txt >> Cmon_allSites_genes_annotations_160818.txt
grep 'Cmoneduloides_001' Cmon_allSites_genes_annotations_160818.txt | awk -F'\t' '{print $16}' | sort -u | wc -l > Cmon_allSites_genes_annotations_numberUniqueGenes_160818.txt

#get genes with warnings about multiple stop codons
awk -F'\t' '$18 ~ "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS" {print $16}' Cmon_allSites_genes_annotations_160818.txt | sort -u > Cmon_allSites_genes_annotations_warningsMultiStop_160818.txt

#remove these genes from main table
grep -vFwf Cmon_allSites_genes_annotations_warningsMultiStop_160818.txt Cmon_allSites_genes_annotations_160818.txt > Cmon_allSites_genes_annotations_noMultiStop_160819.txt
awk -F'\t' '{print $16}' Cmon_allSites_genes_annotations_noMultiStop_160819.txt | sort -u | wc -l Cmon_allSites_genes_annotations_noMultiStop_numberUniqueGenes_160819.txt

#filter the table using custom python scripts
# 1. set genotypes with depth <3 to missing
python /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mk_tests/filter_table_min_depth3_set_missing.py Cmon_allSites_genes_annotations_noMultiStop_160819.txt Cmon_allSites_genes_annotations_noMultiStop_minDP3_160825.txt
awk -F'\t' '{print $16}' Cmon_allSites_genes_annotations_noMultiStop_minDP3_160825.txt | sort -u | wc -l > Cmon_allSites_genes_annotations_noMultiStop_minDP3_numberUniqueGenes_160825.txt

#2. remove lines with >3 genotypes missing and with identical genotypes (i.e. keep only truly variable sites)
python /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mk_tests/filter_table_max_missing2_monomorphic.py  Cmon_allSites_genes_annotations_noMultiStop_minDP3_160825.txt Cmon_allSites_genes_annotations_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt
awk -F'\t' '{print $16}' Cmon_allSites_genes_annotations_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt | sort -u | wc -l > Cmon_allSites_genes_annotations_noMultiStop_minDP3_maxMiss3_onlyVar_numberUniqueGenes_160825.txt


