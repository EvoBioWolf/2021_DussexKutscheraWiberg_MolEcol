### Fisher's exact test (incl. table merging) ###
# originally R version 3.3.2
# failed to install crucial package "FSA" on milou, moved to local system, using R version 3.4.2

mkt <- read.csv("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/mkt_table_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.txt", sep="")

# change the order of the columns
mkt[c(1,4,2,5,3)] -> mkt_sort

# run a Fisher's exact test for each column
get_fisher <- function(df){
  mat <- matrix(as.numeric(df[c(2:5)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(c(df[1], f$p.value))
}

fishers <- apply(mkt_sort, 1,  get_fisher) #apply(X, MARGIN, FUN, ...) where: X is an array or matrix; MARGIN is a variable that determines whether the function is applied over rows (MARGIN=1), columns (MARGIN=2), or both (MARGIN=c(1,2)); FUN is the function to be applied.

# convert the result matrix into a data frame and merge with mkt
fishers_t <- t(fishers) # transpose matrix
fishers_df <- as.data.frame(fishers_t) # convert matrix > data frame
names(fishers_df)[names(fishers_df)=="V2"] <- "Fishers.p.value" # rename the column with the p.values
as.numeric(as.character(fishers_df$Fishers.p.value)) -> fishers_df$Fishers.p.value
merge(mkt_sort, fishers_df, by = "gene") -> mkt_analysed1

# Benjamini-Hochberg correction (FDR) for multiple testing
# order the raw p values
mkt_analysed1 = mkt_analysed1[order(mkt_analysed1$Fishers.p.value),]

# check
library(FSA)
headtail(mkt_analysed1)

# Perform p-value adjustments and add to data frame
mkt_analysed1$BH = p.adjust(mkt_analysed1$Fishers.p.value, method = "BH")


### Neutrality index (NI) (Rand&Kann 1996) ###
#add 1 to all cells if any cell has a count of 0
mkt_sort[apply(mkt_sort[, -1], MARGIN = 1, function(x) any(x == 0)), ] -> mkt_transf
mkt_transf$pN + 1 -> pN
mkt_transf$N.dN + 1 -> N.dN
mkt_transf$pS + 1 -> pS
mkt_transf$S.dS + 1 -> S.dS
data.frame(mkt_transf$gene, pN, N.dN, pS, S.dS) -> mkt_transf
names(mkt_transf)[names(mkt_transf)=="mkt_transf.gene"] <- "gene" # rename the column with gene names

#extract cells where all cells have a count of at least 1 to be merged with transformed data frame
mkt_sort[apply(mkt_sort[, -1], MARGIN = 1, function(x) all(x > 0)), ] -> mkt_untransf
ni_df <- rbind (mkt_transf, mkt_untransf)

#estimate the Neutrality Index NI  
ni_df$NI <- (ni_df$pN / ni_df$N.dN) / (ni_df$pS / ni_df$S.dS)

# merge with mkt_analysed1 (Fisher's exact test)
ni_df_merge <- data.frame(ni_df$gene, ni_df$NI)
names(ni_df_merge)[names(ni_df_merge)=="ni_df.gene"] <- "gene" # rename the column with gene names
names(ni_df_merge)[names(ni_df_merge)=="ni_df.NI"] <- "NI"
merge(mkt_analysed1, ni_df_merge, by = "gene") -> mkt_analysed2
write.table(mkt_analysed2, "/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/MKT_Fishers_NI_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.txt", row.names = FALSE, quote=FALSE)

# Plot Fisher's p values and NI
library("ggplot2")
ggplot(mkt_analysed2, aes(Fishers.p.value)) + geom_histogram(binwidth = 0.01) +
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+ # fontsize
  scale_y_continuous(limits = c(0,20000), breaks = c(0,5000,10000,15000,20000))+
  scale_x_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+
  geom_vline(xintercept=c(0.05), colour="red", linetype="dotted")+
  xlab("p values Fisher's exact test") + guides(fill=FALSE)
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/Fishers_pvalues_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm")

ggplot(mkt_analysed2, aes(Fishers.p.value)) + geom_histogram(binwidth = 0.01) +
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+ # fontsize
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,0.5), breaks = c(0.0,0.1,0.2,0.3,0.4,0.5))+
  geom_vline(xintercept=c(0.05), colour="red", linetype="dotted")+
  xlab("p values Fisher's exact test") + guides(fill=FALSE)
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/Fishers_pvalues_lowFreq_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm")

ggplot(mkt_analysed2, aes(NI)) + geom_histogram(binwidth = 0.1) +
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+ # fontsize
  xlab("NI histogram") + guides(fill=FALSE)+
  geom_vline(xintercept=c(1), colour="red", linetype="dotted")
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/NI_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm") 

ggplot(mkt_analysed2, aes(NI, Fishers.p.value)) + geom_point() + 
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+
  scale_y_continuous(limits = c(0,1.0), breaks = c(0.0,0.2,0.4,0.6,0.8, 1.0))+
  geom_vline(xintercept=c(1), colour="red", linetype="dotted")
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/NI_vs_Fishers_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm")

ggplot(mkt_analysed2, aes(NI)) + geom_histogram(binwidth = 0.1) +
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+ # fontsize
  xlab("NI histogram") + guides(fill=FALSE)+
  scale_x_continuous(limits = c(0,20), breaks = c(0,5,10,15,20))+
  geom_vline(xintercept=c(1), colour="red", linetype="dotted")
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/NI_max20_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm")

ggplot(mkt_analysed2, aes(NI, Fishers.p.value)) + geom_point() +
  theme(panel.background = element_rect(fill = 'white', colour = "black"), text = element_text(size= 17))+
  scale_y_continuous(limits = c(0,1.0), breaks = c(0.0,0.2,0.4,0.6,0.8, 1.0))+
  scale_x_continuous(limits = c(0,20), breaks = c(0,5,10,15,20))+
  geom_vline(xintercept=c(1), colour="red", linetype="dotted")
ggsave("/Users/verenakutschera/Dropbox/Projects/CrowCompGenomicsEBC/Cmoneduloides/MKT/mk_tests/NI_max20_vs_Fishers_pN_pS_noMultiStop_minDP3_maxMiss3_N_S_80_8sp_H_NC_sorted_unique_180118_missingSetTo0.png",width= 11.5, height = 8, units = "cm")

