rm(list=ls())
gc(reset = T,full = T)

library(RColorBrewer)
library(dplyr)
library(tidyr)
library(gplots)
library(Polychrome)
P49 = createPalette(49,  c("#ff0000", "#00ff00", "#0000ff"))
#swatch(P49) for nuclei tRNA isodecoder
P22 = createPalette(22,  c("#ff0000", "#00ff00", "#0000ff"))
#swatch(P22), for MT tRNA isodecoder
###------------------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------------------
#start plot bargraph
dat<-read.delim("Counts/bRMed.combined.cleaned.counts")
protein_sense<-read.delim("Counts/protein.strand")
rownames(protein_sense)<-protein_sense$ID
protein_sense<-protein_sense[c(9:16,1:8),-1]
protein_sense<-t(protein_sense)
protein_sense<-data.frame(prop.table(protein_sense,2))
CDS<-read.delim("Counts/protein_sense.base")
rownames(CDS)<-CDS$ID
CDS<-CDS[c(9:16,1:8),2:5]
CDS<-t(CDS)
CDS<-data.frame(prop.table(CDS,2))

# Profiles
pdf("Figs/profile.pdf",width=11,height=8.5)
par(mfrow=c(3,2))
par(mar=c(3,5,5,2))
# Total RNA profile
scol <- brewer.pal(11, "Set3")[-2]
agg<-aggregate(.~Type3,data=dat[,c(5:21)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
agg.prob$Type=agg$Type
agg.prob<-agg.prob[c(5,6,2,1,4,10,9,3,7,8),]
Type_level<-agg.prob$Type
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=as.vector(agg.prob$Type),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("All RNAs",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

#sncRNA
sncRNA<-dat[dat$Type3=="sncRNA",]
scol <- brewer.pal(9, "Set1")[c(1:3,5:9)]
agg<-aggregate(.~Type,data=sncRNA[,c(3,6:21)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
agg.prob$Type=agg$Type
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=as.vector(agg.prob$Type),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("sncRNAs",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,cex=0.7,y=50,"Reads (%)",srt=90,adj=0.5)

#lncRNA
lncRNA<-dat[dat$Type2=="Antisense"|dat$Type2=="lincRNA"|
            dat$Type2=="Pseudogene"|dat$Type2=="Other lncRNA",]
scol <- brewer.pal(11, "Set3")[c(4:7)]
agg<-aggregate(.~Type3,data=lncRNA[,c(5:21)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
agg.prob$Type=agg$Type
agg.prob<-agg.prob[c(4,1,2,3),]
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=as.vector(agg.prob$Type),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("lncRNAs",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,cex=0.7,y=50,"Reads (%)",srt=90,adj=0.5)

#repeat elements
repeats<-dat[dat$Type3=="Repeats",]
scol <- c("black",brewer.pal(10, "Set3"))
agg<-aggregate(.~Type,data=repeats[,c(3,6:21)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
agg.prob$Type=agg$Type
agg.prob<-agg.prob[c(2,8,5,4,1,3,6,7,9),]
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=as.vector(agg.prob$Type),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Repeat elements",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,cex=0.7,y=50,"Reads (%)",srt=90,adj=0.5)

#protein reads (sense vs anti)
scol <- brewer.pal(11, "Set3")[c(4,5)]
mp<-barplot(as.matrix(protein_sense)*100,col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=c("Sense","Antisense"),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Protein coding gene reads",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)
mtext(rep(paste0("Rep",1:4),4),
      adj=0.5,xpd=TRUE,side=1,line=0.5,cex=0.5,
      at=mp)
mtext(c("SCR-BSA","SCR-PA","KD-BSA","KD-PA"),
      adj=0.5,xpd=TRUE,side=1,line=1.5,cex=0.7,
      at=c(median(mp[1:4]),median(mp[5:8]),median(mp[9:12]),median(mp[13:16])))

#CDS etc
scol <- brewer.pal(5, "Set1")[-3]
mp<-barplot(as.matrix(CDS)*100,col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=c("CDS","UTR","Intron","Intergenic"),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Protein coding gene reads (sense)",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
mtext(rep(paste0("Rep",1:4),4),
      adj=0.5,xpd=TRUE,side=1,line=0.5,cex=0.5,
      at=mp)
mtext(c("SCR-BSA","SCR-PA","KD-BSA","KD-PA"),
      adj=0.5,xpd=TRUE,side=1,line=1.5,cex=0.7,
      at=c(median(mp[1:4]),median(mp[5:8]),median(mp[9:12]),median(mp[13:16])))
par(xpd=T)
text(x=-0.3,y=50,"Bases (%)",srt=90,adj=0.5,cex=0.7,)
dev.off()

pdf("Figs/tRNA_profile.pdf",width=11,height=8.5)
par(mfcol=c(3,2))
par(mar=c(3,5,5,2))
# Nuclear tRNA isoacceptor
tRNA<-dat[dat$Type=="tRNA",]
tmp<-separate(tRNA,col = Name,into = "Iso",sep = "-",extra = "drop")
agg<-aggregate(.~Iso,data=tmp[,c(2,6:21)],FUN=sum)
row.names(agg)<-agg$Iso

agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=P22,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Nuclear tRNAs (Isoacceptor)",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

# Nuclear tRNA isodecoder
agg<-aggregate(.~Name,data=tRNA[,c(2,6:21)],FUN=sum)
row.names(agg)<-agg$Name
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=P49,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.3),
            xlim=c(0,3.5))
mtext("Nuclear tRNAs (Isodecoder)",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

### nuclear tRNA fragments
agg<-aggregate(.~Type2,data=tRNA[,c(4,6:21)],FUN=sum)
row.names(agg)<-agg$Type2
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
scol <- brewer.pal(8, "Set1")[-3]
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Nuclear tRNA fragments",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
mtext(rep(paste0("Rep",1:4),4),
      adj=0.5,xpd=TRUE,side=1,line=0.5,cex=0.5,
      at=mp)
mtext(c("SCR-BSA","SCR-PA","KD-BSA","KD-PA"),
      adj=0.5,xpd=TRUE,side=1,line=1.5,cex=0.7,
      at=c(median(mp[1:4]),median(mp[5:8]),median(mp[9:12]),median(mp[13:16])))
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)


# MT tRNA isoacceptor
tRNA<-dat[dat$Type=="tRNA (MT)",]
tmp<-tRNA
tmp$Name<-substr(tmp$Name,4,5)
agg<-aggregate(.~Name,data=tmp[,c(2,6:21)],FUN=sum)
row.names(agg)<-agg$Name

agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=P22,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Mitochondrial tRNAs (Isoacceptor)",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

# MT tRNA isodecoder
agg<-aggregate(.~Name,data=tRNA[,c(2,6:21)],FUN=sum)
row.names(agg)<-agg$Name
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=P22,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Mitochondrial tRNAs (Isodecoder)",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

### MT tRNA fragments
agg<-aggregate(.~Type2,data=tRNA[,c(4,6:21)],FUN=sum)
row.names(agg)<-agg$Type2
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:17)]),2))
scol <- brewer.pal(8, "Set1")[-3]
mp<-barplot(as.matrix(agg.prob[,c(1:16)]*100),col=scol,axes=F,
            names.arg=rep(NA,16),
            width=0.15,
            space=0.3,
            legend.text=row.names(agg.prob),
            adj=0,args.legend=list(x=3.1,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7),
            xlim=c(0,3.5))
mtext("Mitochondrial tRNA fragments",line=0.25,at=mean(mp),cex=1)
axis(1,labels=NA,at=c(0,3.2),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
axis(1,labels=NA,at=c(mp[1],mp[4]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[5],mp[8]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[9],mp[12]),lwd=1,lwd.ticks=0,line=0.5)
axis(1,labels=NA,at=c(mp[13],mp[16]),lwd=1,lwd.ticks=0,line=0.5)
mtext(rep(paste0("Rep",1:4),4),
      adj=0.5,xpd=TRUE,side=1,line=0.5,cex=0.5,
      at=mp)
mtext(c("SCR-BSA","SCR-PA","KD-BSA","KD-PA"),
      adj=0.5,xpd=TRUE,side=1,line=1.5,cex=0.7,
      at=c(median(mp[1:4]),median(mp[5:8]),median(mp[9:12]),median(mp[13:16])))
par(xpd=T)
text(x=-0.3,y=50,"Reads (%)",srt=90,adj=0.5,cex=0.7)

dev.off()
