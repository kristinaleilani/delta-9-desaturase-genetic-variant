library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggpubr)
library(cowplot)
library(vcfR)
setwd("~/Documents/Ds9")
#------- assembling gene counts table 
counts=read.table('feature_counts_out.tsv',header = T,sep="\t")
genes=counts[ ,1:6]
counts=counts[,-c(1:6)]
bams=colnames(counts)
bams=gsub("allbams/|\\.trim.+","",bams)
head(counts)
row.names(counts)=genes$Geneid
# colnames(counts)=bams


#-------------- retaining genes with mean count>=20

mns=rowMeans(counts)
table(mns>=20)
counts=counts[mns>=20,]


#-------------- setting conditions (origin and transplant site)

sample = colnames(counts)
keep = sample[grep('ub|3m', sample,invert=T)]
counts=counts[,keep]
head(counts)
dim(counts)

sam=sub(".bam","",colnames(counts))
num=sub("[KO]+","",sam)
ori=sub("[KO][0-9]+","",sam)
tra=gsub("^[KO]|[0-9]+$","",sam)
conds=data.frame(cbind(sam,ori,tra))
conds$ind=paste(ori,num,sep="")

dim(counts)
dim(conds)


#---------- Differential gene expression (DESeq2)

dds=DESeqDataSetFromMatrix(counts,
                           colData = conds,
                           design = ~ind+tra) # correcting for genotype

vsd=assay(vst(dds))


dds=DESeq(dds)

tra.s=results(dds,contrast=c("tra","O","K"))
summary(tra.s)
plotMA(tra.s, main="Transplant site- O vs. K", alpha=0.05, colSig = "red")

dim(vsd)

#---------- Volcano plots

vol <- data.frame(gene=row.names(tra.s), pval=-log10(tra.s$padj), lfc=tra.s$log2FoldChange)
# remove na
vol <- na.omit(vol)
# set upper and lower threshold
vol <- mutate(vol, color=case_when(
  vol$lfc > 0 & vol_data2$pval > 1.3 ~ "Increased",
  vol$lfc < 0 & vol_data2$pval > 1.3 ~ "Decreased",
  vol$pval < 1.3 ~ "nonsignificant"))
vol$color[vol_data2$gene == "LOC114946994"] <- "DFAD9"
match_list<-c("LOC114946994")
vol <-vol %>%
  arrange(gene %in% match_list)
vol_plot <- ggplot(vol, aes(x=lfc, y=pval, color=color))
plot_grid(vol_plot + ggtitle(label="Transplant site- O vs. K") +
            geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                               values=c(Increased="#FF66FF", Decreased="#6666FF",
                                        nonsignificant="darkgray", DFAD9="#FF0000")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Transplant site- O vs. K"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), byrow = TRUE, nrow = 2)




#------------- split models (paired by genet)

o2k.c=counts[,conds$ori=="O"]
o2k.meta=conds[conds$ori=="O",]
k2o.c=counts[,conds$ori=="K"]
k2o.meta=conds[conds$ori=="K",]
natives.c=counts[,conds$ori==conds$tra]
natives.meta=conds[conds$ori==conds$tra,]


#-------------  Orpheus to Keppel --------
o2k.dd=DESeqDataSetFromMatrix(o2k.c,
                              colData = o2k.meta,
                              design = ~ind+tra)

o2k.dd=DESeq(o2k.dd)
o2k.vsd=as.data.frame(assay(vst(o2k.dd)))

o2k.ge=results(o2k.dd,contrast=c("tra","K","O")) 
# Comparison will be conducted with Keppel as the numerator, and Orpheus as the denominator
# Genes having positive fold changes will be higher in K transplants when compared to O transplants, and these will appear on the right hand side of the volcano plot.

summary(o2k.ge)

vol_o2k <- data.frame(gene=row.names(o2k.ge), pval=-log10(o2k.ge$padj), lfc=o2k.ge$log2FoldChange)
# remove na
vol_o2k <- na.omit(vol_o2k)
# set upper and lower threshold
vol_o2k <- mutate(vol_o2k, color=case_when(
  vol_o2k$lfc > 0 & vol_o2k$pval > 1.3 ~ "Increased",
  vol_o2k$lfc < 0 & vol_o2k$pval > 1.3 ~ "Decreased",
  vol_o2k$pval < 1.3 ~ "nonsignificant"))
vol_o2k$color[vol_o2k$gene == "LOC114946997"] <- "DFAD9"
match_list<-c("LOC114946997")
vol_o2k <-vol_o2k %>%
  arrange(gene %in% match_list)
vol1 <- ggplot(vol_o2k, aes(x=lfc, y=pval, color=color))

DFAD9.o2k <- plotCounts(o2k.dd, "LOC114946997", "tra", returnData = TRUE)
boxplot(count ~ tra , data=DFAD9.o2k, main = "O to K: Expression of LOC114946997", xlab= "Transplant site")


#-------------  Keppel to Orpheus ------
k2o.dd=DESeqDataSetFromMatrix(k2o.c,
                              colData = k2o.meta,
                              design = ~ind+tra)
k2o.dd=DESeq(k2o.dd)
k2o.vsd=as.data.frame(assay(vst(k2o.dd)))

k2o.ge=results(k2o.dd,contrast=c("tra","K","O"))
summary(k2o.ge)

vol_k2o <- data.frame(gene=row.names(k2o.ge), pval=-log10(k2o.ge$padj), lfc=k2o.ge$log2FoldChange)
# remove na
vol_k2o <- na.omit(vol_k2o)
# set upper and lower threshold
vol_k2o <- mutate(vol_k2o, color=case_when(
  vol_k2o$lfc > 0 & vol_k2o$pval > 1.3 ~ "Increased",
  vol_k2o$lfc < 0 & vol_k2o$pval > 1.3 ~ "Decreased",
  vol_k2o$pval < 1.3 ~ "nonsignificant"))
vol_k2o$color[vol_k2o$gene == "LOC114946997"] <- "DFAD9"
match_list<-c("LOC114946997")
vol_k2o <-vol_k2o %>%
  arrange(gene %in% match_list)
vol2 <- ggplot(vol_k2o, aes(x=lfc, y=pval, color=color))


DFAD9.k2o <- plotCounts(k2o.dd, "LOC114946997", "tra", returnData = TRUE)
boxplot(count ~ tra , data=DFAD9.k2o, main = "K to O: Expression of LOC114946997", xlab= "Transplant site")


#-------------  Natives ------
natives.dd=DESeqDataSetFromMatrix(natives.c,
                                  colData = natives.meta,
                                  design = ~tra)
natives.dd=DESeq(natives.dd)
natives.ge=results(natives.dd,contrast=c("tra","K","O"))
summary(natives.ge)

vol_natives <- data.frame(gene=row.names(natives.ge), pval=-log10(natives.ge$padj), lfc=natives.ge$log2FoldChange)
# remove na
vol_natives <- na.omit(vol_natives)
# set upper and lower threshold
vol_natives <- mutate(vol_natives, color=case_when(
  vol_natives$lfc > 0 & vol_natives$pval > 1.3 ~ "Increased",
  vol_natives$lfc < 0 & vol_natives$pval > 1.3 ~ "Decreased",
  vol_natives$pval < 1.3 ~ "nonsignificant"))
vol_natives$color[vol_natives$gene == "LOC114946997"] <- "DFAD9"
match_list<-c("LOC114946997")
vol_natives <-vol_natives %>%
  arrange(gene %in% match_list)
vol3 <- ggplot(vol_natives, aes(x=lfc, y=pval, color=color))

DFAD9.natives <- plotCounts(natives.dd, "LOC114946997", "tra", returnData = TRUE)
boxplot(count ~ tra , data=DFAD9.natives, main = "Natives: Expression of LOC114946997", xlab= "Transplant site")


library(cowplot)
plot_grid(vol1 +  ggtitle(label="Orpheus to Keppel") +
            geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                               values=c(Increased="#FF66FF", Decreased="#6666FF",
                                        nonsignificant="darkgray", DFAD9="#FF0000")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Orpheus to Keppel"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), 
          vol2 +  ggtitle(label="Keppel to Orpheus") +
            geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                               values=c(Increased="#FF66FF", Decreased="#6666FF",
                                        nonsignificant="darkgray", DFAD9="#FF0000")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Keppel to Orpheus"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), 
          vol3 + ggtitle(label="Natives") +
            geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                               values=c(Increased="#FF66FF", Decreased="#6666FF",
                                        nonsignificant="darkgray", DFAD9="#FF0000")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Natives"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), byrow = TRUE, nrow = 3)

#save(vsd,conds,counts,o2k.ge,k2o.ge,natives.ge,file="Ds9_splitModels.RData")


#-------------  plotting genets over time
o2k.dat<-DFAD9.o2k
all.sam=sub(".bam","",rownames(DFAD9.o2k))
all.num=sub("[KO]+","",all.sam)
o2k.dat=cbind(o2k.dat, all.num)
gg.o2k=ggplot(o2k.dat, aes(tra, count, color = all.num))+
  geom_line(aes(group = all.num))+
  geom_point()+
  ggtitle("Orpheus to Keppel") +
  xlab("Transplant site")+
  labs(color = "Genet")


k2o.dat<-DFAD9.k2o
all.sam=sub(".bam","",rownames(DFAD9.k2o))
all.num=sub("[KO]+","",all.sam)
k2o.dat=cbind(k2o.dat, all.num)
gg.k2o=ggplot(k2o.dat, aes(tra, count, color = all.num))+
  geom_line(aes(group = all.num))+
  geom_point()+ 
  ggtitle("Keppel to Orpheus") +
  xlab("Transplant site")+
  labs(color = "Genet")

library(ggpubr)
ggarrange(gg.o2k, gg.k2o)


#-------------  Combined box plot
all.dat<-cbind(k2o.vsd, o2k.vsd)
all.dat2<-cbind(as.data.frame(t(all.dat)), conds)
all.dat2<-all.dat2[c("tra", "ori", "LOC114946997")]

library(forcats)
ggplot(all.dat2, aes(fct_rev(tra), LOC114946997, fill=fct_rev(ori)))+
  geom_boxplot()+
  scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  labs(fill='Origin site')+
  scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  ggtitle("by origin") +
  ylab("vst transformed counts")+
  xlab("Transplant site")+
  theme_bw()






######---------------- Contrasting by genotype (using significant SNP from 2brad genotypes)
VCF=read.vcfR('GBRall.vcf') 
d2<-as.data.frame(extract.gt(VCF, element = "GT", return.alleles = TRUE))
d2$locus<-row.names(d2)
d2[1:4,1:3] # Look at alleles
d2[d2$locus == "chr9_3997448_T_C", ] # Keppel is mostly C
colnames(d2)=sub(".bam","",colnames(d2))
colnames(d2)=sub("a","",colnames(d2))
d3 = d2[colnames(d2) %in% conds$ind]
d3=as.data.frame(t(d3))
d3 <- tibble::rownames_to_column(d3, "ind")
d4 = d3[,c("ind","chr9_3997448_T_C")] # gts associated with our significant SNP of interest
conds2<-right_join(conds,d4, by="ind")
names(conds2)[5]<-"genotype"
conds2$genotype<-gsub("\\/","", conds2$genotype)
colnames(counts)=sub(".bam","",colnames(counts))
counts2 = counts[colnames(counts) %in% conds2$sam]
conds2 = conds2[!grepl("O8.1", conds2$ind),]
conds2$ori<-as.character(conds2$ori)


#------ Orpheus to Keppel
o2kr.c=counts2[,conds2$ori=="O"]
o2kr.meta=conds2[conds2$ori=="O",]
o2kr.dd=DESeqDataSetFromMatrix(o2kr.c,
                               colData = o2kr.meta,
                               design = ~genotype+tra)
o2kr.dd=DESeq(o2kr.dd)
o2k.vsdt<-as.data.frame(t(o2k.vsd))
o2k.vsdt$sam<-rownames(o2k.vsdt)
o2k.vsdt$sam=sub(".bam","",o2k.vsdt$sam)
o2k.vsd.loci<- o2k.vsdt %>% dplyr::select(c('sam', 'LOC114946997'))
rownames(o2k.vsd.loci)<-o2k.vsd.loci$sam
Ds9.o2kr2 <-merge(o2k.vsd.loci, conds2, by=c("sam"))

#------- Keppel to Orpheus
k2or.c=counts2[,conds2$ori=="K"]
k2or.meta=conds2[conds2$ori=="K",]
k2or.dd=DESeqDataSetFromMatrix(k2or.c,
                               colData = k2or.meta,
                               design = ~genotype+tra)
k2or.dd=DESeq(k2or.dd)
k2o.vsdt<-as.data.frame(t(k2o.vsd))
k2o.vsdt$sam<-rownames(k2o.vsdt)
k2o.vsdt$sam=sub(".bam","",k2o.vsdt$sam)
k2o.vsd.loci<- k2o.vsdt %>% dplyr::select(c('sam', 'LOC114946997'))
Ds9.k2or2 <-merge(k2o.vsd.loci, conds2, by=c("sam"))


allr.dat<-rbind(Ds9.o2kr2, Ds9.k2or2)
allr.dat[allr.dat$ori == "O", "ori"] <- "Orpheus origin" # Renaming for aesthetics
allr.dat[allr.dat$ori == "K", "ori"] <- "Keppel origin"
allr.dat[allr.dat$tra == "O", "tra"] <- "Orpheus Transplant"
allr.dat[allr.dat$tra == "K", "tra"] <- "Keppel Transplant"



CSgt = colorRampPalette(colors = c('skyblue1', 'slateblue1', 'slateblue4')) # create a colorscale for plotting the genotype frequency
Orpheus<-allr.dat %>% dplyr::filter(ori=="Orpheus origin")
a<-ggplot(Orpheus, aes(fct_rev(genotype), LOC114946997, fill=fct_rev(genotype)))+
  geom_boxplot()+
  scale_fill_manual(values=CSgt(3))+
  ylim(5,10)+
  facet_wrap(~fct_rev(tra)) +
  ggtitle("by genotype") +
  ylab("vst transformed counts")+
  xlab("Orpheus origin")+
  theme_bw()

CSgt2 = colorRampPalette(colors = c('slateblue1', 'slateblue4')) # create a colorscale for plotting the genotype frequency
Keppel<-allr.dat %>% dplyr::filter(ori=="Keppel origin")
b<-ggplot(Keppel, aes(fct_rev(genotype), LOC114946997,fill=fct_rev(genotype)))+
  geom_boxplot()+
  scale_fill_manual(values=CSgt2(2))+
  facet_wrap(~fct_rev(tra)) +
  ylab("vst transformed counts")+
  xlab("Keppel origin")+
  theme_bw()

ggarrange(a,b)



