all <- read.delim("/proj/b2013182/nobackup/POPseq/mk_tests/Cmon_allSites_genes_annotations_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt", sep="")

### Subset synonymous and missense variants ###
eff2 <- (all$ANN.EFFECT)
syn_mis2 <- subset(all, eff2 %in% c("missense_variant", "synonymous_variant"))

gene_effect2 <- table(syn_mis2$ANN.FEATUREID, syn_mis2$ANN.EFFECT)
as.data.frame.matrix(gene_effect2) -> gene_effect2_df
gene_effects2_filtered <- gene_effect2_df[c("missense_variant", "synonymous_variant")]
write.table(gene_effects2_filtered, "/proj/b2013182/nobackup/POPseq/mk_tests/mkt_pN_pS_noMultiStop_minDP3_maxMiss3_onlyVar_160825.txt")
