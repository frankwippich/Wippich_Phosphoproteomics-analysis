#####Wippich et al. 2020 
####Phsophoproteomics Data Analysis:

#Loading packages
require(stringr)
require(reshape)
require(dplyr)
require(Biobase)
require(vsn)
require(limma)
require(fdrtool)
require(tidyr)
require(seqinr)
require(ggplot2)
require(easycsv)


#start clean
rm(list=ls())

###set working directory where the MS data, D.m. proeteome and FBgeneidentifier are stored
working_directory<-easycsv::choose_dir()
setwd(working_directory)
###Check if all required files are there
if(file.exists("180419_G0122_MH_FW_trypsin_PSMs.csv")&
   file.exists("Dm_proteome_UP000000803.txt")&
   file.exists("fbgn_NAseq_Uniprot_fb_2018_01.RData"))
{
  

  
  #Loading data and annotating experimental conditions
  conditions<- data.frame("tmt.label"=c("126","127","128","129","130","131"), "rep"=c("rep1","rep2","rep1","rep2","rep3","rep3"),"condition"=c("fed","fed","stvd","stvd","stvd","fed"), "File.name"=rep("180419_G0122_MH_FW_trypsin_PSMs",6))
  
  data<-NULL
  for(file in unique(conditions$File.name)){
    x<-read.csv(file.path(paste0(file,".csv")))
    x<-subset(x,!grepl("###",Protein.Group.Accessions))
    x$File.name<-file
    data<-rbind(data,x)
    rm(x)
  }
  
  #Call Modifications
  data$phosphosite <- (grepl("Phospho",data$Modifications))
  data$phosphosite[data$phosphosite =="TRUE"] <- "P"
  data$phosphosite[data$phosphosite =="FALSE"] <- "noP"
  data$peptide_entry_ID <- paste(rownames(data),data$Protein.Group.Accessions,data$Sequence,data$phosphosite, sep = "_")
  data$peptide_entry_ID_top3 <- paste(data$Protein.Group.Accessions,data$Sequence,data$phosphosite, sep = "_")
  rm(file)
  
  
  
  
  ########################
  ##########Proteome changes:
  ###############################################
  
  #Building an expression set
  ##Restructuring data
  data2<-data[,c("peptide_entry_ID","peptide_entry_ID_top3","File.name","Modifications","Sequence","phosphosite",
                 "Spectrum.File",
                 "IonScore",grep("signal_sum",names(data),value = T))]
  #######################
  # FILTER
  data2<-subset(data2,IonScore >15)
  data2$IonScore<-NULL
  
  
  # only Flow Through of the TiO2 is used to calc proteme change
  data2<-subset(data2,grepl("_FT_",Spectrum.File))
  # only non phosphorylted peptides of the FT are used to calc proteme change
  Data_proteome<- subset(data2, phosphosite=="noP")
  Data_proteome$phosphosite<-NULL
  
  #use major Uniprot ID for LIMMA 
  Data_proteome$UniprotID<- str_split_fixed(Data_proteome$peptide_entry_ID, "_", 3)[,2]
  Data_proteome$UniprotID.group<-Data_proteome$UniprotID
  Data_proteome$UniprotID<- str_split_fixed(Data_proteome$UniprotID, "-", 2)[,1]
  Data_proteome$UniprotID<-str_split_fixed(Data_proteome$UniprotID, ";", 2)[,1]
  
  Data_proteome$UniprotID[Data_proteome$UniprotID==""]<-NA
  
  
  # mTop3
  Data_proteome2<-Data_proteome[,c("UniprotID","signal_sum_126")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_126)
  Data_proteome2<-Data_proteome2 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2.1<-Data_proteome[,c("UniprotID","signal_sum_127")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_127)
  Data_proteome2.1<-Data_proteome2.1 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2<-merge(Data_proteome2,Data_proteome2.1)
  rm(Data_proteome2.1)
  
  Data_proteome2.1<-Data_proteome[,c("UniprotID","signal_sum_128")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_128)
  Data_proteome2.1<-Data_proteome2.1 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2<-merge(Data_proteome2,Data_proteome2.1)
  rm(Data_proteome2.1)
  
  Data_proteome2.1<-Data_proteome[,c("UniprotID","signal_sum_129")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_129)
  Data_proteome2.1<-Data_proteome2.1 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2<-merge(Data_proteome2,Data_proteome2.1)
  rm(Data_proteome2.1)
  
  Data_proteome2.1<-Data_proteome[,c("UniprotID","signal_sum_130")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_130)
  Data_proteome2.1<-Data_proteome2.1 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2<-merge(Data_proteome2,Data_proteome2.1)
  rm(Data_proteome2.1)
  
  Data_proteome2.1<-Data_proteome[,c("UniprotID","signal_sum_131")] %>% group_by(UniprotID) %>% top_n(n=3,signal_sum_131)
  Data_proteome2.1<-Data_proteome2.1 %>% group_by(UniprotID) %>% summarise_each(funs(median))
  
  Data_proteome2<-merge(Data_proteome2,Data_proteome2.1)
  rm(Data_proteome2.1)
  
  
  #merge Filename
  Data_proteome2$File.name<-conditions$File.name[1]
  
  
  mdata.P<-reshape::melt(Data_proteome2,id.vars = c("UniprotID","File.name"),variable_name = c("tmt.label"))
  
  mdata.P$tmt.label<-gsub("([a-z,A-Z]+_)+","",mdata.P$tmt.label)
  mdata.P$tmt.label<-as.character(mdata.P$tmt.label)
  mdata.P$measurement<-gsub("signal_sum","",mdata.P$tmt.label)
  mdata.P<-merge(mdata.P,conditions)
  
  cdata.P<-cast(mdata.P,formula=UniprotID~measurement+condition+rep+tmt.label,value = "value",fun.aggregate = mean,na.rm=T)
  cdata.P$found.in.samples<- rep("3",length(cdata.P$UniprotID))
  
  cdata.P<-data.frame(cdata.P)
  
  ##Constructing assay data
  raw_data.P<-cdata.P[,c("UniprotID",grep("_rep",names(cdata.P),value=T))]
  rownames(raw_data.P)<-with(raw_data.P,paste0(UniprotID))
  raw_data.P$UniprotID<-NULL
  # raw_data$description<-NULL
  raw_data.P<-as.data.frame(raw_data.P)
  
  ##Constructing metadata
  metadata.P<-data.frame(col.name=names(raw_data.P))
  metadata.P$tmt.label<-str_split_fixed(metadata.P$col.name, "_", 4)[,4]
  metadata.P$rep<-str_split_fixed(metadata.P$col.name, "_", 4)[,3]
  metadata.P<-merge(metadata.P,conditions,sort=F)
  metadata.P$ID<-metadata.P$col.name
  
  ##Constructing feature data
  feat_data.P<-cdata.P[,c("UniprotID","found.in.samples")]
  rownames(feat_data.P)<-feat_data.P$UniprotID
  rownames(metadata.P)<-metadata.P$ID
  colnames(raw_data.P)<-metadata.P$ID
  
  ##Transformation of raw signal_sum data
  # The log2 is computed from the raw signal intensities. Furthermore, Infinite values are transformed into missing (NA) values.
  raw_data_m.P<-log2(as.matrix(raw_data.P))
  raw_data_m.P[is.infinite((raw_data_m.P))]<-NA
  raw_data_m.P[is.na((raw_data_m.P))]<-0
  
  ##Constructing the Expression set
  raw_dataE.P <- ExpressionSet(assayData = raw_data_m.P, 
                               featureData = AnnotatedDataFrame(feat_data.P),
                               phenoData = AnnotatedDataFrame(metadata.P))
  
  #Normalization
  #The vsn package from Wolfgang Huber is used to apply a variance stabilization normalization method on the log2 raw data.
  vsn.fit.P<-vsn2(2^exprs(raw_dataE.P))
  norm_dataE.P<-raw_dataE.P
  #meanSdPlot((vsn.fit.P))
  exprs(norm_dataE.P)<-predict(vsn.fit.P,2^exprs(norm_dataE.P))
  rm(vsn.fit.P)
  
  #LIMMA analysis to identify differential proteins
  
  ##Construction of design matrix and defining conditions that will be compared
  condition<-factor(pData(norm_dataE.P)$condition)
  
  condition_d<-model.matrix(~0+condition)
  colnames(condition_d)<-gsub("condition","",colnames(condition_d))
  to_label<-c("stvd - fed")
  
  ##Fitting the model
  limma_comparison<-eBayes(contrasts.fit(lmFit(norm_dataE.P,design=condition_d),
                                         makeContrasts(contrasts = to_label,levels=condition_d)))
  rm(condition_d,condition)
  
  
  
  ##Identification of top hits
  data2_proteome<-NULL
  for(comp in to_label){
    res <- topTable(limma_comparison, sort.by = "t",  coef = comp, number = Inf)
    fdr_res <- fdrtool(res$t, plot =F, verbose = F)
    res$local.p.value<-fdr_res$pval
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    data2_proteome<-rbind(data2_proteome,res)
  }
  
  data2_proteome$fdr<-p.adjust(data2_proteome$local.p.value,method = "fdr")
  rm(res,limma_comparison,comp,to_label)
  
  p.value.limit=0.05
  fc.limit<-2
  data2_proteome$hit<-with(data2_proteome,ifelse((logFC>log2(fc.limit)|(logFC<log2(fc.limit)))&P.Value<p.value.limit,T,F))
  
  
  ### Now bring back ID groups 
  data2_proteome$logFC.Proteome<-data2_proteome$logFC
  data2_proteome$P.Value.Proteome<-data2_proteome$P.Value
  data2_proteome<-unique(merge(data2_proteome,unique(Data_proteome[,c("UniprotID","UniprotID.group")]),all.x = TRUE))
  
  
  #identify FLybase Gene ID
  if(!exists("FBgn_conversion")){
    load("fbgn_NAseq_Uniprot_fb_2018_01.RData")
    names(FBgn_conversion)[names(FBgn_conversion)=="gene_symbol"] <- "gene_name"
    names(FBgn_conversion)[names(FBgn_conversion)=="primary_FBgn."] <- "Flybase_gene_ID"
    names(FBgn_conversion)[names(FBgn_conversion)=="nucleotide_accession"] <- "UniprotID"
    FBgn_conversion<-FBgn_conversion[,c("gene_name","Flybase_gene_ID","UniprotID")]
  }
  
  data2_proteome<-merge(data2_proteome,FBgn_conversion, by="UniprotID")
  
  write.csv(data2_proteome,paste("data_proteome",
                                 # Sys.time(),
                                 ".csv",sep = ""))
  
  
  
  ###############################################
  ##Restructuring data for Phospho-site analysis
  ###############################################
  data2<-data[,c("peptide_entry_ID","File.name","Modifications","Sequence","phosphosite",
                 "Spectrum.File",
                 "IonScore",grep("signal_sum",names(data),value = T))]
  
  # FILTER
  data2<-subset(data2,IonScore >15 & phosphosite == "P")
  data2$IonScore<-NULL
  # only TiO2 bound
  data2<-subset(data2,!grepl("_FT_",Spectrum.File))
  
  ###########Put back Gene_ID
  data2$UniprotID<- str_split_fixed(data2$peptide_entry_ID, "_", 3)[,2]
  data2$UniprotID<- str_split_fixed(data2$UniprotID, "-", 2)[,1]
  data2$unique <- !grepl(";",data2$UniprotID) 
  data2.not.unique<-subset(data2, unique == "FALSE")
  data2.not.unique.multiplied<-separate_rows(data2.not.unique,UniprotID,sep=";",convert = T)
  data2<-rbind(subset(data2, unique == "TRUE"), data2.not.unique.multiplied)
  
  # new unique Peptide_ID
  data2$peptide_entry_ID<-paste(seq(1:length(data2$peptide_entry_ID)),data2$UniprotID, data2$Sequence,sep="_")
  
  #######
  if(!exists("FBgn_conversion")){
    load("fbgn_NAseq_Uniprot_fb_2018_01.RData")
    names(FBgn_conversion)[names(FBgn_conversion)=="gene_symbol"] <- "gene_name"
    names(FBgn_conversion)[names(FBgn_conversion)=="primary_FBgn."] <- "Flybase_gene_ID"
    names(FBgn_conversion)[names(FBgn_conversion)=="nucleotide_accession"] <- "UniprotID"
    FBgn_conversion<-FBgn_conversion[,c("gene_name","Flybase_gene_ID","UniprotID")]
  }
  
  data2<-merge(data2,FBgn_conversion, by="UniprotID")
  
  
  #Load and compare Proteoms from FASTA/TXT
  if(!exists("FASTA")){
    fastaFile <- read.fasta("Dm_proteome_UP000000803.txt",as.string = TRUE,forceDNAtolower = FALSE)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    FASTA <- data.frame(seq_name, sequence)
    FASTA$UniprotID<- str_split_fixed(FASTA$seq_name, "_", 3)[,2]
    rm(fastaFile,sequence,seq_name)
  }
  
  
  # map peptides to protein
  FASTA<-FASTA %>% dplyr::rename(sequence.FASTA = sequence)
  data2<-merge(data2,FASTA[,c("UniprotID","sequence.FASTA")])
  
  data2<-data2%>% group_by(peptide_entry_ID) %>% mutate(pep_start= regexpr(toupper(Sequence),sequence.FASTA))
  data2$pep_start[data2$pep_start== -1]<-NA
  
  data2$pep_stop <- data2$pep_start+nchar(as.character(data2$Sequence))-1
  data2$prot_length <- nchar(as.character(data2$sequence.FASTA))
  
  
  
  #Identify the Phospho-position
  data2<-separate_rows(data2,Modifications,sep=";",convert = T) ###Expand -> subset P
  data2<-subset(data2, grepl("Phospho",data2$Modifications))    ###-> subset P
  data2$Modifications<-gsub(" ", "", data2$Modifications) 
  data2$Phospho_position_Nr<-as.numeric(gsub("[A-Z,a-z,()]","",data2$Modifications))   
  data2$Phospho_position_Nr<-mapply('+',data2$Phospho_position_Nr,(data2$pep_start-1))
  data2$Phospho_position<-paste(substr(as.character(data2$Modifications),1,1),data2$Phospho_position_Nr,sep="")
  
  
  #####unique IDs
  data2$ID <- paste(data2$UniprotID,data2$Phospho_position,sep="_")
  
  ######Merge all Sites by median top3
  data3<-data2[,c("ID",grep("signal_sum",names(data),value = T))]
  data3$ID<-as.factor(data3$ID)
  data4<-data.frame(unique(data3[,"ID"]))
  names(data4)[1] <- "ID"
  
  # mTop3
  data4<-data3[,c("ID","signal_sum_126")] %>% group_by(ID) %>% top_n(n=3,signal_sum_126)
  data4<-data4 %>% group_by(ID) %>% summarise_each(funs(median))
  data4.1<-data3[,c("ID","signal_sum_127")] %>% group_by(ID) %>% top_n(n=3,signal_sum_127)
  data4.1<-data4.1 %>% group_by(ID) %>% summarise_each(funs(median))
  data4<-merge(data4,data4.1)
  rm(data4.1)
  data4.1<-data3[,c("ID","signal_sum_128")] %>% group_by(ID) %>% top_n(n=3,signal_sum_128)
  data4.1<-data4.1 %>% group_by(ID) %>% summarise_each(funs(median))
  data4<-merge(data4,data4.1)
  rm(data4.1)
  data4.1<-data3[,c("ID","signal_sum_129")] %>% group_by(ID) %>% top_n(n=3,signal_sum_129)
  data4.1<-data4.1 %>% group_by(ID) %>% summarise_each(funs(median))
  data4<-merge(data4,data4.1)
  rm(data4.1)
  data4.1<-data3[,c("ID","signal_sum_130")] %>% group_by(ID) %>% top_n(n=3,signal_sum_130)
  data4.1<-data4.1 %>% group_by(ID) %>% summarise_each(funs(median))
  data4<-merge(data4,data4.1)
  rm(data4.1)
  data4.1<-data3[,c("ID","signal_sum_131")] %>% group_by(ID) %>% top_n(n=3,signal_sum_131)
  data4.1<-data4.1 %>% group_by(ID) %>% summarise_each(funs(median))
  data4<-merge(data4,data4.1)
  rm(data4.1)
  
  data4<-unique(merge(data4,data2[,c("ID","File.name")],by="ID"))
  
  ###########
  mdata<-reshape::melt(data4,id.vars = c("ID","File.name"),variable_name = c("tmt.label"))
  
  mdata$tmt.label<-gsub("([a-z,A-Z]+_)+","",mdata$tmt.label)
  mdata$tmt.label<-as.character(mdata$tmt.label)
  mdata$measurement<-gsub("signal_sum","",mdata$tmt.label)
  mdata<-merge(mdata,conditions)  
  
  cdata<-cast(mdata,formula=ID~measurement+condition+rep+tmt.label,value = "value",fun.aggregate = mean,na.rm=T)
  cdata$found.in.samples<- rep("3",length(cdata$ID))
  cdata<-data.frame(cdata)
  
  
  ##Constructing assay data
  raw_data<-cdata[,c("ID",grep("_rep",names(cdata),value=T))]
  rownames(raw_data)<-with(raw_data,paste0(ID))
  raw_data$ID<-NULL
  raw_data<-as.data.frame(raw_data)
  
  ##Constructing metadata
  metadata<-data.frame(col.name=names(raw_data))
  metadata$tmt.label<-str_split_fixed(metadata$col.name, "_", 4)[,4]
  metadata$rep<-str_split_fixed(metadata$col.name, "_", 4)[,3]
  metadata<-merge(metadata,conditions,sort=F)
  metadata$ID<-metadata$col.name
  
  ##Constructing feature data
  feat_data<-cdata[,c("ID","found.in.samples")]
  rownames(feat_data)<-feat_data$ID
  rownames(metadata)<-metadata$ID
  colnames(raw_data)<-metadata$ID
  
  ##Transformation of raw signal_sum data
  # The log2 is computed from the raw signal intensities. Furthermore, Infinite values are transformed into missing (NA) values.
  raw_data_m<-log2(as.matrix(raw_data))
  raw_data_m[is.infinite((raw_data_m))]<-NA
  raw_data_m[is.na((raw_data_m))]<-0
  
  ##Constructing the Expression set
  raw_dataE <- ExpressionSet(assayData = raw_data_m, 
                             featureData = AnnotatedDataFrame(feat_data),
                             phenoData = AnnotatedDataFrame(metadata))
  
  #Normalization
  #The vsn package from Wolfgang Huber is used to apply a variance stabilization normalization method on the log2 raw data.
  vsn.fit<-vsn2(2^exprs(raw_dataE))
  norm_dataE<-raw_dataE
  #meanSdPlot((vsn.fit))
  exprs(norm_dataE)<-predict(vsn.fit,2^exprs(norm_dataE))
  rm(vsn.fit)
  
  #LIMMA analysis to identify differential proteins
  
  ##Construction of design matrix and defining conditions that will be compared
  condition<-factor(pData(norm_dataE)$condition)
  
  condition_d<-model.matrix(~0+condition)
  colnames(condition_d)<-gsub("condition","",colnames(condition_d))
  to_label<-c("stvd - fed")
  
  ##Fitting the model
  limma_comparison<-eBayes(contrasts.fit(lmFit(norm_dataE,design=condition_d),
                                         makeContrasts(contrasts = to_label,levels=condition_d)))
  rm(condition_d,condition)
  
  ##Identification of top hits
  limma_res<-NULL
  # limma_old<-NULL
  for(comp in to_label){
    res <- topTable(limma_comparison, sort.by = "t",  coef = comp, number = Inf)
    fdr_res <- fdrtool(res$t, plot =F, verbose = F)
    res$local.p.value<-fdr_res$pval
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    limma_res<-rbind(limma_res,res)
  }
  
  limma_res$fdr<-p.adjust(limma_res$local.p.value,method = "fdr")
  rm(res,limma_comparison,comp,to_label)
  
  p.value.limit=0.05
  fc.limit<-2
  limma_res$hit<-with(limma_res,ifelse((logFC>log2(fc.limit)|(logFC< -log2(fc.limit)))&P.Value<p.value.limit,T,F))
  
  
  ###########Put back UniprotID,other info
  
  # limma_res$UniprotID<- str_split_fixed(limma_res$peptide_entry_ID, "[0-9]_", 2)[,1]
  limma_res$UniprotID<- str_split_fixed(limma_res$ID, "_", 2)[,1]
  limma_res$Phospho_position<- str_split_fixed(limma_res$ID, "_", 2)[,2]
  limma_res$Phospho_position_Nr<-gsub("[A-Z]","",limma_res$Phospho_position)
  
  info<-unique(data2[,c("UniprotID","prot_length")])
  limma_res<-merge(limma_res,info)
  rm(info)
  
  if(!exists("FBgn_conversion")){
    load("fbgn_NAseq_Uniprot_fb_2018_01.RData")
    names(FBgn_conversion)[names(FBgn_conversion)=="gene_symbol"] <- "gene_name"
    names(FBgn_conversion)[names(FBgn_conversion)=="primary_FBgn."] <- "Flybase_gene_ID"
    names(FBgn_conversion)[names(FBgn_conversion)=="nucleotide_accession"] <- "UniprotID"
    FBgn_conversion<-FBgn_conversion[,c("gene_name","Flybase_gene_ID","UniprotID")]
  }
  
  limma_res<-merge(limma_res,FBgn_conversion, by="UniprotID",all.x = TRUE)
  
  
  #########
  #Merge with Proteome Analysis
  data3_proteome <- separate_rows(data2_proteome,UniprotID.group,sep=";",convert = T)
  limma_res<-merge(limma_res,data3_proteome[,c("UniprotID","logFC.Proteome","P.Value.Proteome")], all.x=TRUE)
  limma_res<-unique(limma_res)
  
  
  ################
  # normalize the P-sites with Proteome changes
  
  # when proteome data available:
  limma_res.protnorm<-subset(limma_res,!is.na(limma_res$logFC.Proteome)) # when proteome data available:
  # when logFC.Proteome<0:
  limma_res.protnorm.lesszero<-subset(limma_res.protnorm,logFC.Proteome<0)
  limma_res.protnorm.lesszero$norm_logFC<-limma_res.protnorm.lesszero$logFC+(-1*limma_res.protnorm.lesszero$logFC.Proteome)
  # when logFC.Proteome>0:
  limma_res.protnorm.greatzero<-subset(limma_res.protnorm,logFC.Proteome>0)
  limma_res.protnorm.greatzero$norm_logFC<-limma_res.protnorm.greatzero$logFC-(limma_res.protnorm.greatzero$logFC.Proteome)
  # when proteome data is not available:
  limma_res.notprotnorm<-subset(limma_res,is.na(limma_res$logFC.Proteome)) 
  limma_res.notprotnorm$norm_logFC<-limma_res.notprotnorm$logFC
  #merge em back
  limma_res<-rbind(limma_res.protnorm.greatzero,limma_res.protnorm.lesszero,limma_res.notprotnorm)  
  
  rm(limma_res.notprotnorm,limma_res.protnorm,limma_res.protnorm.greatzero,limma_res.protnorm.lesszero)
  
  #call hit's again
  limma_res$hit<-(limma_res$norm_logFC>log2(fc.limit)|limma_res$norm_logFC< -log2(fc.limit))&limma_res$P.Value<p.value.limit
  limma_res$hit.signif<-with(limma_res,ifelse(P.Value<0.05 & (norm_logFC < -1 |norm_logFC >1) ,"*",""))
  
  ##########get sequence of Phospho-site back
  data.seq<-data2[,c("ID","UniprotID","peptide_entry_ID","Modifications","Sequence","Phospho_position")]
  data.seq<- merge(data.seq,FASTA,by ="UniprotID", all.x = TRUE)
  data.seq$Phospho_position_Nr<-as.numeric(gsub("[A-Z]","",data.seq$Phospho_position)) 
  data.seq$Sequence_around_P<- NULL
  data.seq$Sequence_around_P<- substr(as.character(data.seq$sequence),
                                      data.seq$Phospho_position_Nr-5,    #start
                                      data.seq$Phospho_position_Nr+4)    #stop
  
  limma_res$Sequence_around_P.x<-NULL
  limma_res$Sequence_around_P.y<-NULL
  limma_res$Sequence_around_P<-NULL
  
  limma_res<- merge(limma_res,data.seq[,c("ID","Sequence_around_P")],by ="ID", all.x = TRUE)
  limma_res<-unique(limma_res)
  
  
  
  # remove duplicates due to duplicated UniprotID
  limma_res$ID2<-paste(limma_res$gene_name,limma_res$Sequence_around_P,sep = "_")
  limma_res$ID3<-paste(limma_res$gene_name,limma_res$Sequence_around_P, limma_res$norm_logFC,sep = "_")
  
  limma_res$duplicate<-duplicated(limma_res$ID3)
  limma_res_with_dupl<-limma_res
  limma_res=limma_res_with_dupl[!duplicated(limma_res_with_dupl$ID3, fromLast = TRUE), ]
  
  ##################################################################################################
  ###Kinase Consensus Site Mapping
  ##################################################################################################
  {
    
    CK2_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][D][A-Z][ED][A-Z]")
    
    limma_res$CK2_consens <- grepl(paste(CK2_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    CK1_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][A-Z][A-Z][ST][A-Z]")
    
    limma_res$CK1_consens <- grepl(paste(CK1_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    Cdc2_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][P][A-Z][RK][A-Z]")
    
    limma_res$Cdc2_consens <- grepl(paste(Cdc2_consens,collapse = "|"),
                                    limma_res$Sequence_around_P )
    
    PKA_consens<- c("[A-Z][A-Z][A-Z][R][A-Z][ST][A-Z][A-Z][A-Z][A-Z]" ,"[A-Z][A-Z][R][RK][A-Z][ST][A-Z][A-Z][A-Z][A-Z]","[A-Z][KR][A-Z][A-Z][ST][A-Z][A-Z][A-Z][A-Z]")
    
    limma_res$PKA_consens <- grepl(paste(PKA_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    AKT_consens<- c("[R][A-Z][R][A-Z][A-Z][ST][A-Z][A-Z][A-Z][A-Z]")
    
    limma_res$AKT_consens <- grepl(paste(AKT_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    GSK3_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][A-Z][A-Z][ST]")
    
    limma_res$GSK3_consens <- grepl(paste(GSK3_consens,collapse = "|"),
                                    limma_res$Sequence_around_P )
    
    Aurora_consens<- c("[A-Z][A-Z][A-Z][RK][A-Z][ST][LVI][A-Z][A-Z]")
    
    limma_res$Aurora_consens <- grepl(paste(Aurora_consens,collapse = "|"),
                                      limma_res$Sequence_around_P )
    
    CDK1_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][P][A-Z][KR]","[A-Z][A-Z][A-Z][A-Z][A-Z][ST][P][KR][A-Z]")
    
    limma_res$CDK1_consens <- grepl(paste(CDK1_consens,collapse = "|"),
                                    limma_res$Sequence_around_P )
    
    ERKandMAPK_consens<- c("[A-Z][A-Z][A-Z][PV][A-Z][ST][P][A-Z][A-Z][A-Z]")
    
    limma_res$ERKandMAPK_consens <- grepl(paste(ERKandMAPK_consens,collapse = "|"),
                                          limma_res$Sequence_around_P )
    
    AMPK_consens<- c("[L][A-Z][R][A-Z][A-Z][ST][A-Z][A-Z][A-Z][A-Z]")
    
    limma_res$AMPK_consens <- grepl(paste(AMPK_consens,collapse = "|"),
                                    limma_res$Sequence_around_P )
    
    TOR_consens<- c("[A-Z][A-Z][A-Z][A-Z][A-Z][ST][FPLWYV][A-Z][A-Z][A-Z]")
    
    limma_res$TOR_consens <- grepl(paste(TOR_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    S6K_consens<- c("[RK][A-Z][R][A-Z][A-Z][ST][A-Z][A-Z][A-Z][A-Z]")
    
    limma_res$S6K_consens <- grepl(paste(S6K_consens,collapse = "|"),
                                   limma_res$Sequence_around_P )
    
    
    CAMK2_consens<- c("[A-Z][A-Z][R][A-Z][A-Z][ST][A-Z][A-Z][A-Z][A-Z]","[A-Z][A-Z][R][A-Z][A-Z][ST][V][A-Z][A-Z][A-Z]")
    
    limma_res$CAMK2_consens <- grepl(paste(CAMK2_consens,collapse = "|"),
                                     limma_res$Sequence_around_P )
  }
  
  
  ##################################################################################################
  ### Insulin/TOR pathway analysis Fig. 6b
  ##################################################################################################
  
  InR.list<-c("RpS6",
              "Akt1",
              "chico",
              "foxo",
              "gig",
              "Pdk1",
              "Pi3K68D",
              "S6kII",
              "Thor",
              "eIF4B",
              "eIF4G1",
              "eIF4G2",
              "eIF4H1",
              "Tor")
  
  limma_res$InR.pathway <- data.matrix(limma_res$gene_name %in% InR.list)
  
  TEST_InR<-t.test(limma_res$norm_logFC[limma_res$InR.pathway==FALSE],limma_res$norm_logFC[limma_res$InR.pathway==TRUE],paired=FALSE)

  
  plot_InR<-ggplot(data=limma_res,aes(InR.pathway,norm_logFC,fill=InR.pathway))+
    ggtitle("Fig. 6b: InR pathway components")+
    geom_jitter(aes(color=InR.pathway.handpicked), size=0.5)+
    geom_boxplot(size=0.2, position=pd,outlier.shape = NA,alpha=0.2)+
    scale_fill_manual(values = c("grey","#F0B949"))+
    scale_color_manual(values = c("grey","#F0B949"))+
    labs(x=paste("\n","p = ",signif(TEST_hand$p.value, digits=3)))+
    theme_classic()
  
  ##################################################################################################
  ### HDAC1 phosphorylation sites Fig. 7a
  ##################################################################################################
  i="Q94517"
  
  sub<-subset(limma_res_with_dupl[,c("UniprotID","gene_name","Phospho_position","logFC","P.Value","prot_length")],UniprotID==i)#grepl(i,gene_name))
  sub$data<-"phospho"
  sub$Phospho_position_Nr<-as.numeric(gsub("[A-Z]","",sub$Phospho_position))   
  sub.P<-unique(subset(limma_res_with_dupl[,c("UniprotID","gene_name","logFC.Proteome","P.Value.Proteome")], grepl(i,limma_res_with_dupl$UniprotID)))#grepl(i,gene_name))
  sub.P<-dplyr::rename(sub.P, logFC = logFC.Proteome)
  sub.P<-dplyr::rename(sub.P, P.Value = P.Value.Proteome)
  sub.P$Phospho_position<-"Protein"
  sub.P$Phospho_position_Nr<- -5
  sub.P$prot_length<-sub$prot_length[1]
  sub.P$data<-"protein"
  subbind<-rbind(sub,sub.P)
  subbind$norm_logFC<-subbind$logFC+(-1*subbind$logFC[subbind$data=="protein"])
  subbind$significant[subbind$data=="phospho"]<-subbind$P.Value[subbind$data=="phospho"]<0.05 & (subbind$norm_logFC[subbind$data=="phospho"] < -1 |subbind$norm_logFC[subbind$data=="phospho"] >1) 
  subbind$significant[subbind$data=="protein"]<-subbind$P.Value[subbind$data=="protein"]<0.05 & (subbind$logFC[subbind$data=="protein"] < -1 |subbind$logFC[subbind$data=="protein"] >1) 
  
  
  RoundUp <- function(from,to) ceiling(from/to)*to
  
  HDAC1_plot<- ggplot(subbind)+
    geom_segment(data=subset(subbind, data=="protein"),aes(x=0,xend=prot_length,y=norm_logFC,yend=norm_logFC),size=2,color="lightblue")+
    geom_segment(aes(x=1, xend=374,y=0,yend=0),size=2)+
    geom_text(aes(x=150, y=0,label="Histone dacetylase domain"),vjust = -0.8,hjust = 0.2 ,size=4)+
    geom_point(data=subset(subbind, data=="phospho"),aes(Phospho_position_Nr,norm_logFC,fill=significant),size=2,shape=21)+
    geom_text(data=subset(subbind, data=="phospho"), aes(Phospho_position_Nr,norm_logFC,label=Phospho_position),hjust=-0.15, vjust=-0.15,size = 2.5, color="firebrick")+
    xlim(-5,RoundUp(subbind$prot_length[1], 100) )+
    scale_fill_manual(values=c("#ADCDE3","#E41A1C"))+
    ylim(-3.5,
         3.5)+
    xlab("peptide position within protein") +
    ggtitle("Fig. 7a",paste(unique(limma_res_with_dupl$gene_name[limma_res_with_dupl$UniprotID==i]),i,unique(limma_res_with_dupl$Flybase_gene_ID[limma_res_with_dupl$UniprotID==i]),sep="_"))+
    theme_bw()
  
  
  
  #################
  #Write Figures to  PDF
  pdf(file= paste("Figures.pdf"), width = 6, height = 5,paper = "A4")
  plot_InR
  HDAC1_plot
  dev.off()
  
  
  #################
  #Write Results
  results<-limma_res[,c("UniprotID",
                        "Flybase_gene_ID",
                        "gene_name",
                        "Phospho_position",
                        "prot_length",
                        "hit",
                        "norm_logFC",
                        "P.Value",
                        "logFC.Proteome",
                        "P.Value.Proteome",
                        "Sequence_around_P",
                        "InR.pathway",
                        grep("consens",names(limma_res),value=T))]
  write.csv(results,"results.csv")
  
}else{print("ERROR: Missing Files")}
