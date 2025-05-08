rm(list=ls())
gc(reset = T,full = T)

library(dplyr)
library(DESeq2)
library(VennDiagram)
library(grid)
library(RColorBrewer)
flog.threshold(ERROR)
library(ggplot2)
library(tidyr)
library(sva)
library(tidyverse)
dat<-read.delim("Counts/bRMed.combined.cleaned.counts")
### use batch removed counts
## ribosomal RNA removed
dat<-dat[dat$Type3!="rRNA",]
dat<-dat[dat$Type!="rRNA (MT)",]

coldata<-data.frame("Cell"=rep(c("SCR","KD"),each=8),
                    "Treatment"=rep(rep(c("BSA","PA"),each=4),2),
                    "Sample"=colnames(dat)[6:21])
coldata$Comb<-paste0(coldata$Cell,"_",coldata$Treatment)
coldata$Cell<-factor(coldata$Cell,levels = c("SCR","KD"))
coldata$Treatment<-factor(coldata$Treatment,levels = c("BSA","PA"))
coldata$Comb<-factor(coldata$Comb,levels = c("SCR_BSA","SCR_PA","KD_BSA","KD_PA"))
### make DESEq2 object
'''
P_dds <- DESeqDataSetFromMatrix(countData = dat[,c(6:21)],
                                     colData = coldata,design = ~Cell+Treatment+Cell:Treatment)
P_dds <- DESeq(P_dds,parallel = T,test = "LRT",reduced = ~1)
saveRDS(P_dds,"DEseqObj/JS_batchRMed_LRT.rds")
rld <- rlog(P_dds, blind=TRUE)
saveRDS(rld,"DEseqObj/dds_LRT_rld_bRMed.rds")
'''

### LRT
dds_lrt<-readRDS("DEseqObj/JS_batchRMed_LRT.rds")
res_LRT <- results(dds_lrt)

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="ID") %>% 
  as_tibble() %>% 
  filter(padj <= 0.05 )


# Subset results for faster cluster finding (for classroom demo purposes)
clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) 
# Get sig gene lists
sigLRT_genes <- clustering_sig_genes %>% 
  pull(ID)

rld<-readRDS("DEseqObj/dds_LRT_rld_bRMed.rds")
rld_mat <- assay(rld) 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

rld_mat <- data.frame(t(apply(rld_mat, 1, cal_z_score)))
rld_mat<-cbind(rld_mat,dat[,1:5])

bcol<- brewer.pal(10, "Set3")[c(-2,-5)]
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)
top1000<-rld_mat[clustering_sig_genes$ID,]
top1000$Type3<-factor(top1000$Type3,levels=c("Protein coding","Pseudogene","lincRNA","Other lncRNA","tRNA","sncRNA","MT","Repeats"))
top1000$Name[top1000$Type=="miRNA"]<-paste(top1000$Name[top1000$Type=="miRNA"],top1000$Type2[top1000$Type=="miRNA"],sep=":")
top1000$Name[top1000$Type=="tRNA"]<-paste(top1000$Name[top1000$Type=="tRNA"],top1000$Type2[top1000$Type=="tRNA"],sep=":")
pdf("Figs/LRT_heatmap.pdf",height=15,width=8)
hm<-heatmap.2(as.matrix(top1000[,1:16]),labRow = top1000$Name,dendrogram = "row",
          scale="none",margins = c(10, 20),breaks=seq(-3,3,0.06),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,symkey=F,Colv=F,colRow = bcol[top1000$Type3],
          col =heatPalette,RowSideColors=bcol[top1000$Type3],
          cexRow = 0.5, cexCol=0.75,
          lhei=c(1,16), lwid=c(1,10), keysize=0.5, key.par = list(cex=0.5),
          ColSideColors=c(rep("gray",4),rep("orange",4),rep("skyblue",4),
                          rep("gold",4)))
legend("topright",
       legend = levels(top1000$Type3),
       col = bcol, 
       lty= 1,             
       lwd = 5,           
       cex=.5)
dev.off()


hc <- as.hclust( hm$rowDendrogram)
group<-cutree( hc, h=6.3 )

## 4 groups
###
hc<-merge(top1000,data.frame(group),by=0)
hc<-hc[,c(2:18,23)]
hc$SCR_BSA<-rowMeans(hc[,1:4])
hc$SCR_PA<-rowMeans(hc[,5:8])
hc$KD_BSA<-rowMeans(hc[,9:12])
hc$KD_PA<-rowMeans(hc[,13:16])
hc<-hc[,17:22]
hc <- hc %>%
  pivot_longer(
    cols = `SCR_BSA`:`KD_PA`, 
    names_to = "type",
    values_to = "value",
    values_transform = ~ as.numeric(gsub(",", "", .x))
  ) %>% as.data.frame
hc$type<-factor(hc$type,levels = c("SCR_BSA","SCR_PA","KD_BSA","KD_PA"))
hc$group<-factor(hc$group)
hc$ID<-factor(hc$ID)
pdf("Figs/deg_pattern.pdf")
ggplot(hc, aes(x = type, y = value, group=type, color = type)) + 
  geom_boxplot(alpha = 0, outlier.size = 0, outlier.shape = NA) + 
  geom_point(alpha = 0.4, size = 1, position = position_jitterdodge(dodge.width = 0.9)) + 
  stat_smooth(aes(x = type, y = value, group=group,color=type), se = FALSE, method = "lm", 
              formula = y ~ poly(x, 3)) +
  geom_line(aes(group = ID,color=type), 
            alpha = 0.1) +
  facet_wrap(~group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust = 0.5)) + ylab("Z-score of gene abundance") + xlab("")
dev.off()

