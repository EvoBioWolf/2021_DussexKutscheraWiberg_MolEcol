#!/bin/bash

# edit the output from R
sed -i.bak 's/"//g' mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt
sed -i.bak 's/Cmoneduloides_001-//g' mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt

# prepare the polymorphism data table for the McDonald-Kreitman test
grep -v '^missense_variant' mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt | sort -k1 > mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_sorted_160826.txt

# add header line "gene pN pS"
echo -e "gene pN pS" | cat - mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_sorted_160826.txt >  mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_sorted_header_160826.txt
