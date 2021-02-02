###
# Analysis of the crow codeml results
#
#
###
# Clear all
rm(list = ls(all=TRUE))

#
# Load libraries #####
library(ggplot2)
library(stringr)
source("~/RData/RScripts/ggplot_theme.R")
library(reshape2)
library(dplyr)
library(plyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
# ####

#
species5<-c("Csple","Cfru","Cmon","3sp_Ccorone","Cdau","Tgut")
# Load data ####
nspecies<-"5species"
seqs<-"_mskVars"
A_data<- read.table(paste("~/Projects/crows/paml_results/",
                          nspecies,
                          "_A_results_pamlinfiles",seqs,".txt",sep=""),
                    sep=",",header = FALSE)
A_data$V1<-gsub(" ","",A_data$V1)
head(A_data)
# Add header to dataset
if(nspecies=="8species"){
  # 8 species dataset: 3sp_Ccornix,Ccorx,Chaw,Cfru,Csple,Ctas,Cmon,Cdau,Tgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species8,"_dN",sep=""),
                      paste(species8,"_dS",sep=""),
                      paste(species8,"_w",sep=""),
                      "lnL","np","model")
  # Load the 8species complete genes
  complg<-read.table(paste("crows_paml/",
                   nspecies,
                   "_complete_genes_pamlinfiles.list",sep=""),
             sep=",",header = FALSE)
  colnames(complg)<-c("gene")
  # Subset
  A_data<-A_data[A_data$gene %in% complg$gene,]
}else if(nspecies=="7species"){
  # 7 species dataset: 3sp_Ccornix,Cdau,Cfru,Csple,Cmon,Ctas,CcorxTgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species7,"_dN",sep=""),
                      paste(species7,"_dS",sep=""),
                      paste(species7,"_w",sep=""),
                      "lnL","np","model")
}else if(nspecies=="5species"){
  # 5 species dataset: 3sp_Ccornix,Cdau,Cfru,Csple,Cmon,Tgut ###
  colnames(A_data)<-c("gene","stop","N","S",
                      paste(species5,"_dN",sep=""),
                      paste(species5,"_dS",sep=""),
                      paste(species5,"_w",sep=""),
                      "lnL","np","model")
}
head(A_data)
# ####
Null_data<- read.table(paste("~/Projects/crows/paml_results/",
                             nspecies,
                             "_Null_results_pamlinfiles",seqs,".txt",sep=""),
                       sep=",",header = FALSE)
Null_data$V1<-gsub(" ","",Null_data$V1)
# Add header to dataset
if(nspecies=="5species"){
  # 5species dataset:
  colnames(Null_data)<-c("gene","stop","N","S","dN",
                         "dS","w","lnL","np","model",
                         paste(species5,"_dN",sep=""),
                         paste(species5,"_dS",sep=""))
  }else if(nspecies == "7species"){
    # 7species dataset:
    colnames(Null_data)<-c("gene","stop","N","S","dN",
                           "dS","w","lnL","np","model",
                           paste(species7,"_dN",sep=""),
                           paste(species7,"_dS",sep=""))
  }else if(nspecies == "8species"){
    # 8species dataset:
    colnames(Null_data)<-c("gene","stop","N","S","dN",
                           "dS","w","lnL","np","model",
                           paste(species8,"_dN",sep=""),
                           paste(species8,"_dS",sep=""))
    # Subset
    Null_data<-Null_data[Null_data$gene %in% complg$gene,]
  }

# ####
nrow(A_data)
nrow(Null_data)

head(A_data)
head(Null_data)
#
# Distribution of omega
# How many with w > 10
nrow(Null_data[Null_data$w > 10,])
nrow(A_data[A_data$Cmon_w > 10,])
nrow(A_data[A_data$`3sp_Ccorone_w` > 10,])

ggplot()+
  geom_histogram(data=Null_data[Null_data$w < 10,],aes(w))

A_data_m<-melt(data=A_data,measure.vars = paste(species5,"_w",sep=""))
head(A_data_m)
ggplot()+
  geom_histogram(data=A_data_m[A_data_m$value < 10,],aes(value))+
  facet_grid(variable~.)

# Perform LRT ####
A_data$p <- vector(length=nrow(A_data))
for (g in A_data$gene){
  np1 <- A_data$np[A_data$gene == g]
  np2 <- Null_data$np[Null_data$gene == g]
  lnL1 <- A_data$lnL[A_data$gene == g]
  lnL2 <- Null_data$lnL[Null_data$gene == g]
  LR <- 2*(lnL1 - lnL2)
  df <- np1-np2
  A_data$p[A_data$gene == g]<-pchisq(LR,df,lower.tail = FALSE)
  Null_data$p[Null_data$gene == g]<-pchisq(LR,df,lower.tail = FALSE)
}
# ####

#
# Subset from A_data those with significant LRT ####
# Convert p-values to q-values
qobj <- qvalue(A_data$p)
plot(qobj)
A_data$qvalues <- qobj$qvalues
Null_data$qvalues <- qobj$qvalues
head(A_data)

# Which are the "significant" genes (q-value < 0.05)
sig_A_data <- A_data[A_data$qvalues < 1,]
# Remove the spaces in the gene names
sig_A_data$gene <- gsub(" ","",sig_A_data$gene)
# How many
nrow(sig_A_data)
sig_A_data[,1]
# Show the data
head(sig_A_data)
tail(sig_A_data)
# nrow(sig_A_data[sig_A_data$p < bonf,])
sig_A_data
A_data[A_data$gene=="DTNBP1 ",]

beak_genes<-c("POU1F1","IGF2R","BMP4","")

head(A_data[A_data$p < 0.05,])
sig_A_data<-A_data[A_data$p < 0.05,]
sig_A_data[which(sig_A_data$gene %in% genes),]

# Write significant results to a table
write.csv(sig_A_data,
          paste("~/Projects/crows/paml_results/",
                nspecies,
                seqs,
                "_sig_A_data.tab",
                sep=""),
          row.names = FALSE,quote = FALSE)

# ####


# Compare data to Nic's results
# Load Nic's Null data:
nic_null_data<-read.table("~/Projects/crows/nic_paml/5sp_Ns_Null_results_pamlinfiles.txt",
                          header=FALSE,sep=",")
colnames(nic_null_data)<-c("nic_gene","nic_stop","nic_N","nic_S","nic_dN","nic_dS","nic_w","nic_lnl","nic_np","nic_model")
nic_null_data$nic_gene<-gsub(".fa","",nic_null_data$nic_gene)
nic_null_data$nic_gene<-gsub(" ","",nic_null_data$nic_gene)
head(nic_null_data)
nrow(nic_null_data)


# Load Nic's A data:
nic_A_data<-read.table("~/Projects/crows/nic_paml/5sp_Ns_A_results_pamlinfiles.txt",
                          header=FALSE,sep=",")
nic_5sp<-c("Cmon","Cfru","Csple","Ccorone","Cdau","Tgut")
colnames(nic_A_data)<-c("nic_gene","nic_stop","N","S",
                        paste(nic_5sp,"_dN",sep=""),
                        paste(nic_5sp,"_dS",sep=""),
                        paste(nic_5sp,"_w",sep=""),
                        "nic_lnl","nic_np","nic_model")
nic_A_data$nic_gene<-gsub(".fa","",nic_A_data$nic_gene)
nic_A_data$nic_gene<-gsub(" ","",nic_A_data$nic_gene)

head(nic_A_data)

# Load the ortholog table and add Tgut gene names to the Null dataset
nic_orths<-read.table("~/Projects/crows/nic_paml/orthologs_genes_Tgutt_noties.abc.txt",
                      header=FALSE,sep="\t")
colnames(nic_orths)<-c("Tgut","Cmon","score")
nic_orths$Tgut<-gsub("Tgut\\|","",nic_orths$Tgut)
nic_orths$Cmon<-gsub("Cmone\\|Cmoneduloides_001-","",nic_orths$Cmon)
head(nic_orths)
# How many names
length(nic_orths$Cmon)
length(nic_orths$Tgut)
# How many unique names
length(unique(nic_orths$Cmon))
length(unique(nic_orths$Tgut))

nic_null_data$tgut_gene<-vector(length=nrow(nic_null_data))
nic_A_data$tgut_gene<-vector(length=nrow(nic_A_data))
for(name in nic_null_data$nic_gene){
  if(name %in% nic_orths$Cmon){
    tgut_name<-nic_orths$Tgut[nic_orths$Cmon==name]
    nic_null_data$tgut_gene[nic_null_data$nic_gene==name]<-tgut_name
    nic_A_data$tgut_gene[nic_A_data$nic_gene==name]<-tgut_name
  }else{
    nic_null_data$tgut_gene[nic_null_data$nic_gene==name]<-NA
    nic_A_data$tgut_gene[nic_A_data$nic_gene==name]<-NA
  }
}

# Subset to the genes that have Tgut ID
nic_null_data_filt<-nic_null_data[!(is.na(nic_null_data$tgut_gene)),]
nic_A_data_filt<-nic_null_data[!(is.na(nic_A_data$tgut_gene)),]
head(nic_null_data_filt,n=50)
head(nic_A_data_filt,n=50)
nrow(nic_null_data_filt)
nrow(nic_A_data_filt)


# Perform LRT ####
nic_A_data$p <- vector(length=nrow(nic_A_data))
nic_null_data$p <- vector(length=nrow(nic_null_data))
for (g in nic_A_data$nic_gene){
  np1 <- nic_A_data$nic_np[nic_A_data$nic_gene == g]
  np2 <- nic_null_data$nic_np[nic_null_data$nic_gene == g]
  lnL1 <- nic_A_data$nic_lnl[nic_A_data$nic_gene == g]
  lnL2 <- nic_null_data$nic_lnl[nic_null_data$nic_gene == g]
  LR <- 2*(lnL1 - lnL2)
  df <- np1-np2
  nic_A_data$p[nic_A_data$nic_gene == g]<-pchisq(LR,df,lower.tail = FALSE)
  nic_null_data$p[nic_null_data$nic_gene == g]<-pchisq(LR,df,lower.tail = FALSE)
}
# ####

#
# Subset from A_data those with significant LRT ####
# Convert p-values to q-values
qobj <- qvalue(nic_A_data$p)
plot(qobj)
nic_A_data$qvalues <- qobj$qvalues
nic_null_data$qvalues <- qobj$qvalues


nic_sig_A_data <- nic_A_data[nic_A_data$qvalues < 1,]

# Match the datasets
length(which(Null_data$gene %in% nic_null_data_filt$tgut_gene))
Null_data$nic_w<-nic_null_data$nic_w[match(Null_data$gene,nic_null_data$tgut_gene)]
Null_data$nic_dn<-nic_null_data$nic_dN[match(Null_data$gene,nic_null_data$tgut_gene)]
Null_data$nic_ds<-nic_null_data$nic_dS[match(Null_data$gene,nic_null_data$tgut_gene)]
Null_data$nic_q<-nic_null_data$qvalues[match(Null_data$gene,nic_null_data$tgut_gene)]
Null_data$nic_p<-nic_null_data$p[match(Null_data$gene,nic_null_data$tgut_gene)]
head(Null_data)
nrow(Null_data[complete.cases(Null_data),])
Null_data_filt<-Null_data[complete.cases(Null_data),]

nic_null_data[nic_null_data$nic_w > 10,]

head(Null_data_filt)
# Spearman Rank Correlation test 
cor.test(Null_data_filt$w,Null_data_filt$nic_w,method = "spearman")

ggplot()+
  geom_point(data=Null_data_filt[Null_data_filt$nic_w < 10,],aes(x=w,y=nic_w))+
  xlab(expression(paste(omega,"(",italic("C. c. cornix")," as Reference)")))+
  ylab(expression(paste(omega,"(",italic("C. moneduloides")," as Reference)")))+
  geom_label(aes(x=0.2,y=2,label="Spearman Rank Correlation:\nrho = 0.83\np < 0.001"))+
  my.theme  

# Spearman Rank Correlation test 
cor.test(Null_data_filt$dN,Null_data_filt$nic_dn,method = "spearman")

ggplot()+
  geom_point(data=Null_data_filt[Null_data_filt$nic_w < 10,],aes(x=dN,y=nic_dn))+
  xlab(expression(paste("dN (",italic("C. c. cornix")," as Reference)")))+
  ylab(expression(paste("dN (",italic("C. moneduloides")," as Reference)")))+
  geom_label(aes(x=0.2,y=2,label="Spearman Rank Correlation:\nrho = 0.83\np < 0.001"))+
  my.theme  

# Spearman Rank Correlation test
cor.test(Null_data_filt$dS,Null_data_filt$nic_ds,method = "spearman")

ggplot()+
  geom_point(data=Null_data_filt[Null_data_filt$nic_w < 10,],aes(x=dS,y=nic_ds))+
  xlab(expression(paste("dN (",italic("C. c. cornix")," as Reference)")))+
  ylab(expression(paste("dN (",italic("C. moneduloides")," as Reference)")))+
  geom_label(aes(x=0.2,y=2,label="Spearman Rank Correlation:\nrho = 0.83\np < 0.001"))+
  my.theme  

# Spearman Rank Correlation test 
cor.test(Null_data_filt$qvalues,Null_data_filt$nic_q,method = "spearman")
ggplot()+
  geom_point(data=Null_data_filt[Null_data_filt$nic_w < 10,],aes(x=-log10(qvalues),y=-log10(nic_q)))+
  xlab(expression(paste("-log10(q-values) (",italic("C. c. cornix")," as Reference)")))+
  ylab(expression(paste("-log10(q-values) (",italic("C. moneduloides")," as Reference)")))+
  geom_label(aes(x=30,y=30,label="Spearman Rank Correlation:\nrho = 0.35\np < 0.001"))+
  my.theme  

# Spearman Rank Correlation test 
cor.test(Null_data_filt$p,Null_data_filt$nic_p,method = "spearman")
ggplot()+
  geom_point(data=Null_data_filt[Null_data_filt$nic_w < 10,],aes(x=-log10(p),y=-log10(nic_p)))+
  xlab(expression(paste("-log10(p-values) (",italic("C. c. cornix")," as Reference)")))+
  ylab(expression(paste("-log10(p-values) (",italic("C. moneduloides")," as Reference)")))+
  geom_label(aes(x=30,y=30,label="Spearman Rank Correlation:\nrho = 0.55\np < 0.001"))+
  my.theme  


sig_A_data$nic_Cmon_w<-nic_sig_A_data$Cmon_w[match(sig_A_data$gene,nic_sig_A_data$tgut_gene)]
sig_A_data$nic_3sp_Ccorone_w<-nic_sig_A_data$Ccorone_w[match(sig_A_data$gene,nic_sig_A_data$tgut_gene)]
sig_A_data[complete.cases(sig_A_data),]
sig_A_data
head(sig_A_data)

# Spearman Rank Correlation test 
cor.test(sig_A_data$Cmon_w,sig_A_data$nic_Cmon_w,method = "spearman")
ggplot()+
  geom_point(data=sig_A_data,aes(x=Cmon_w,y=nic_Cmon_w))+
  xlab("Cmon dN/dS (C. c. cornix as Reference)")+
  ylab("Cmon dN/dS (C. moneduloides as Reference")

ggplot()+
  geom_point(data=sig_A_data,aes(x=`3sp_Ccorone_w`,y=nic_3sp_Ccorone_w))+
  xlab("BG dN/dS (C. c. cornix as Reference)")+
  ylab("BG dN/dS (C. moneduloides as Reference")





