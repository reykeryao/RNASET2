rm(list=ls())
gc(reset = T,full = T)

library(dplyr)
library(tidyr)
### ---JS11248 2-fraction dataset counts cleanup-----
### tRFs counts, and  mature/hairpin miRNA counts
### no ERCC counts
counts_loc<-"Counts/"
dat<-read.table(paste0(counts_loc,"combined.counts"),header = T)
### remove tRNA, MT_tRNA, miRNA
dat<-dat[dat$Type!="tRNA",]
dat<-dat[dat$Type!="Mt_tRNA",]
dat<-dat[dat$Type!="miRNA",]
### --- cleanup biotype
dat<-dat[dat$Type!="SP",]
dat<-dat[dat$Type!="ERCC",]
dat$Type[dat$Type=="protein_coding"]<-"Protein coding"
dat$Type[dat$Type=="TR"]<-"Protein coding"
dat$Type[dat$Type=="IG"]<-"Protein coding"
dat$Type[dat$Type=="pseudogene"]<-"Pseudogene"
dat$Type[dat$Type=="miscRNA"]<-"misc RNA"
dat$Type[dat$Type=="lncRNA"]<-"Other lncRNA"
dat$Type[dat$Type=="VTRNA"]<-"VT RNA"
dat$Type[dat$Type=="YRNA"]<-"Y RNA"
dat$Type[dat$Type=="18S_rRNA"]<-"rRNA"
dat$Type[dat$Type=="28S_rRNA"]<-"rRNA"
dat$Type[dat$Type=="5.8S_rRNA"]<-"rRNA"
dat$Type[dat$Type=="5S_rRNA"]<-"rRNA"
dat$Type[dat$Type=="scaRNA"]<-"snoRNA"
dat$Type[dat$Type=="Mt_protein_coding"]<-"Protein coding (MT)"
dat$Type[dat$Type=="Mt_rRNA"]<-"rRNA (MT)"
dat$Type2<-dat$Type
dat$Type2[dat$Type2=="Protein coding (MT)"]<-"MT"
dat$Type2[dat$Type2=="rRNA (MT)"]<-"MT"
dat$Type2[dat$Type2=="7SK"]<-"sncRNA"
dat$Type2[dat$Type2=="7SL"]<-"sncRNA"
dat$Type2[dat$Type2=="snoRNA"]<-"sncRNA"
dat$Type2[dat$Type2=="snRNA"]<-"sncRNA"
dat$Type2[dat$Type2=="VT RNA"]<-"sncRNA"
dat$Type2[dat$Type2=="Y RNA"]<-"sncRNA"
dat$Type2[dat$Type2=="misc RNA"]<-"sncRNA"
dat$Type3<-dat$Type2
dat<-dat[,c(1:3,20:21,12:19,4:11)]
### replace MT tRNA counts
MT_tRNA<-read.table(paste0(counts_loc,"MT_tRNA_frag.counts"),header = T)
MT_tRNA$Name<-MT_tRNA$ID
MT_tRNA$Type<-"tRNA (MT)"
MT_tRNA$Type3<-"MT"
colnames(MT_tRNA)[2]<-"Type2"
MT_tRNA$ID<-paste(MT_tRNA$Name,MT_tRNA$Type2,sep=":")
MT_tRNA<-MT_tRNA[,c(1,19:20,2,21,11:18,3:10)]
dat<-rbind(dat,MT_tRNA)
rm(MT_tRNA)
### replace tRNA counts
tRNA<-read.table(paste0(counts_loc,"tRNA_frag.counts"),header = T)
tRNA$Name<-tRNA$ID
tRNA$Type<-"tRNA"
tRNA$Type3<-"tRNA"
colnames(tRNA)[2]<-"Type2"
tRNA$ID<-paste(tRNA$Name,tRNA$Type2,sep=":")
tRNA<-tRNA[,c(1,19:20,2,21,11:18,3:10)]
dat<-rbind(dat,tRNA)
rm(tRNA)
### replace miRNA
miRNA<-read.table(paste0(counts_loc,"combined_miRNAs.counts"),header = T)
miRNA$Type2<-miRNA$Type
miRNA<-separate(miRNA,col="ID",into = c("ID","Name","Type"),sep=":")
miRNA$Type3<-"sncRNA"
miRNA<-miRNA[,c(1:3,20:21,12:19,4:11)]
miRNA$ID<-paste(miRNA$Name,miRNA$Type2,sep=":")
dat<-rbind(dat,miRNA)
rm(miRNA)
### add repeat elements
Rep<-read.table(paste0(counts_loc,"combined_repeats.counts"),header = T)
colnames(Rep)[1:3]<-c("Name","Type","Type2")
### sub family
### non-LTR retrotransposons
###  LINE--already subtyped
Rep$Type2[Rep$Type=="LINE" & Rep$Type2!="L1" & Rep$Type2!="L2"]<-"Other LINE"
Rep$Type2[Rep$Type=="LINE" & Rep$Type2=="L1" ]<-"Other L1"
Rep$Type2[grep("^L1HS",Rep$Name)]<-"L1HS"
Rep$Type2[grep("^L1PA",Rep$Name)]<-"L1PA2-17"
###  SINE
Rep$Type2[Rep$Type=="SINE"]<-"SINE"
Rep$Type2[grep("AluY|^Alu$",Rep$Name)]<-"AluY"
Rep$Type2[grep("AluS",Rep$Name)]<-"AluS"
Rep$Type2[grep("AluJ",Rep$Name)]<-"AluJ"
Rep$Type2[grep("FAM|FLAM|FRAM",Rep$Name)]<-"Monomeric Alu"
Rep$Type2[grep("MIR|LFSINE|MamSINE1",Rep$Name)]<-"SINE2|tRNA"
Rep$Type2[Rep$Type2=="SINE"]<-"SINE3|5S rRNA"
### composite
Rep$Type2[Rep$Type=="Retroposon"]<-"SVA"

### LTR retrotransposon
###  LTR
Rep$Type2[Rep$Type=="LTR"]<-"LTR"
Rep$Type2[grep("MamGyp|LTR85|LTR88|LTR104",Rep$Name)]<-"Gypsy"
Rep$Type2[grep("^MLT|^HERVL-|^ERVL|^ERVL-|^MST|^THE1",Rep$Name)]<-"ERV3"
Rep$Type2[grep("^HERV16|^LTR16|^LTR47|^LTR18|^LTR66|^HERVL66|^MER68|^LTR40|^HERVL40|^HERVL18",Rep$Name)]<-"ERV3"
Rep$Type2[grep("^HERV74|^MER74|^LTR52|^LTR57|^MER70|^MER76|^LTR53",Rep$Name)]<-"ERV3"
name_tmp<-c("LTR32", "LTR33", "LTR33A", "LTR33B", "LTR33C", "LTR41", "LTR41B", 
            "LTR41C", "LTR42", "LTR50", "LTR55", "LTR62", 
            "LTR67B", "LTR69", "LTR75", "LTR75B", "LTR79", "LTR80A", "LTR80B", 
            "LTR82A", "LTR82B", "LTR83", "LTR84a", "LTR84b", "LTR86A1", 
            "LTR86A2", "LTR86B1", "LTR86B2", "LTR86C", "LTR87", "LTR89", 
            "LTR91", "LTR108d_Mam", "LTR108e_Mam", "MER54", "MER54A", "MER54B", 
            "MER73", "MER74", "MER74A", "MER74B", "MER77", "RMER10B")
Rep$Type2[Rep$Name%in%name_tmp]<-"ERV3"
Rep$Type2[grep("^HERVK",Rep$Name)]<-"ERV2"
name_tmp<-c("LTR14A", "LTR14B","LTR5", "LTR5A","MER9a1", "MER9a2", "MER9a3",
            "LTR13", "LTR13A","LTR22A", "LTR22B", "LTR22B1", "LTR22B2", 
            "LTR22C", "LTR22C0", "LTR22C2","LTR3", "LTR3A", "LTR3B",
            "MER11D","MER11A", "MER11B", "MER11C","LTR14C","LTR14","LTR5B", 
            "LTR5_Hs", "LTR22", "LTR22E", "MER9", "MER9B", "RLTR10B", "RLTR10C")
Rep$Type2[Rep$Name%in%name_tmp]<-"ERV2"
Rep$Type2[Rep$Type2=="LTR"]<-"ERV1"
### DNA transposon
### DNA
Rep$Type2[Rep$Type=="DNA"]<-"DNA"
name_tmp<-c("Eulor5A", "Eulor5B", "Eulor6A", "Eulor6B", "Eulor6C", "Eulor6D", 
            "Eulor6E")
Rep$Type2[Rep$Name%in%name_tmp]<-"Crypton"
Rep$Type2[grep("^Rick",Rep$Name)]<-"MuDR"
Rep$Type2[grep("^UCON29",Rep$Name)]<-"Kolobok"
Rep$Type[grep("^Helitron",Rep$Name)]<-"DNA"
Rep$Type2[grep("^Helitron",Rep$Name)]<-"Helitron"
Rep$Type2[grep("^Merlin",Rep$Name)]<-"Merlin"
name_tmp<-c("LOOPER","Looper", "MER75", "MER75A", "MER75B", "MER85")
Rep$Type2[Rep$Name%in%name_tmp]<-"piggyBac"
name_tmp<-c("EutTc1-N1", "GOLEM","Tigger3", "GOLEM_A","Tigger3a", "GOLEM_B", 
            "GOLEM_C", "HSMAR1", "HSMAR2", "HSTC2", "Kanga1", "Kanga1d", 
            "KANGA2_A","Kanga2_a", "Kanga11a", "MADE1", "MARE1", "MARE10", 
            "MARNA", "MER2", "MER2B", "MER6", "MER6A", "MER6B", "MER6C", 
            "MER8", "MER28","Tigger2a", "MER44A", "MER44B", "MER44C", 
            "MER44D", "MER46C", "MER47B", "MER47C", "MER53", "MER82", 
            "MER104", "MER104A","Kanga1a", "MER104B","Kanga1b", "MER104C",
            "Kanga1c", "MER116", "MER127", "MER132", "MERX", "MamRep137", 
            "MamRep434", "TIGGER1","Tigger1", "TIGGER2","Tigger2", "Tigger3b", 
            "Tigger3c", "Tigger3d", "Tigger4a", "TIGGER5","Tigger5", "TIGGER5A",
            "MER47A", "TIGGER5_A", "TIGGER5_B","Tigger5b", "TIGGER6A",
            "Tigger6a", "TIGGER6B","Tigger6b", "TIGGER7", "Tigger7", "TIGGER8",
            "Tigger8", "TIGGER9", "Tigger9b", "Tigger10", "Tigger12", "Tigger12A",
            "Tigger13a", "Tigger14a", "Tigger15a", "Tigger16a", "Tigger16b",
            "Tigger2b_Pri", "UCON39", "UCON42", "UCON104", "X1_DNA", "X6a_DNA", 
            "X6b_DNA", "X10a_DNA", "X10b_DNA", "X13_DNA", "X25_DNA", "X26_DNA", 
            "X32_DNA", "X33a_DNA", "ZOMBI","Tigger4", "ZOMBI_A", "ZOMBI_B", "ZOMBI_C")
Rep$Type2[Rep$Name%in%name_tmp]<-"Mariner|Tc1"
Rep$Type2[Rep$Type2=="DNA"]<-"hAT"
### simplify type
Rep$Type[Rep$Type=="Low_complexity"]<-"Low complexity"
Rep$Type[Rep$Type=="Simple_repeat"]<-"Simple repeat"
Rep$Type[Rep$Type=="Unknown"]<-"Unclassified"
### simplify type2
Rep$Type2[Rep$Type=="Low complexity"]<-"Low complexity"
Rep$Type2[Rep$Type=="Simple repeat"]<-"Simple repeat"
Rep$Type2[Rep$Type=="Satellite"]<-"Satellite"
Rep$Type2[Rep$Type=="Unclassified"]<-"Unclassified"
Rep$ID<-paste(Rep$Name,Rep$Type,Rep$Type2,sep=":")
### aggregate by ID name
Rep<-aggregate(.~ID,data=Rep[,c(4:20)],FUN = sum)
Rep<-separate(Rep,col="ID",into=c("Name","Type","Type2"),sep=":",remove=F)
Rep$Type3="Repeats"
Rep<-Rep[,c(1:4,21,13:20,5:12)]
dat<-rbind(dat,Rep)
rownames(dat)<-dat$ID
write.table(dat,paste0(counts_loc,"combined.cleaned.counts"),sep="\t",row.names=T,quote=F)
rm(Rep,name_tmp)


### take-care of batch effect
library(sva)
batch<-t(as.vector(rep(c(1,1,2,2),4)))
bRmCounts<-ComBat_seq(counts =as.matrix(dat[,c(6:21)]),batch = batch)
bRmCounts<-cbind(dat[,1:5],bRmCounts)
write.table(bRmCounts,paste0(counts_loc,"bRMed.combined.cleaned.counts"),quote=F,sep="\t",row.names=T)
