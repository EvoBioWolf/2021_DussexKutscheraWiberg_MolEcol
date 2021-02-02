### Table preparation ###
poly <- read.csv("/proj/b2013182/nobackup/POPseq/mk_tests/mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_sorted_header_160826.txt", sep="")
fixed <- read.csv("/proj/b2013182/nobackup/POPseq/mk_tests/mkt_table_N_S_80_7sp_sorted_unique.txt", sep="")

# merge the two tables by "gene"
merge(x = fixed, y = poly, by = "gene", all = TRUE) -> mkt

# check genes with missing data
subset(mkt, is.na(N.dN)) -> N.dN.NA
subset(mkt, is.na(pN)) -> pN.NA
write.table(N.dN.NA, "/proj/b2013182/nobackup/POPseq/mk_tests/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_N_S.txt", row.names = FALSE)
write.table(pN.NA, "/proj/b2013182/nobackup/POPseq/mk_tests/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_missing_genes_pN_pS.txt", row.names = FALSE)

# export merged table
write.table(mkt, "/proj/b2013182/nobackup/POPseq/mk_tests/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_7sp_sorted_unique_180118.txt", row.names = FALSE)

