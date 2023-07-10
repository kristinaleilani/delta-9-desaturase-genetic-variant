setwd("~/Documents/Ds9") 
library(vcfR)
library(dplyr)
library(ggplot2)
library(qqman)
library(stringr)



###------ Import output from Bayescan
# Based on scripts from https://rpubs.com/lbenestan/outlier by Laura Benestan
bayescan=read.table("GBRa.baye_fst.txt") 
SNPb=read.table("bsa.sites",header=FALSE)
bayescan=cbind(SNPb, bayescan) 
bayescan <- bayescan %>%
  tidyr::unite("SNP", V1:V2)
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

# Change the value of the Q_VALUE column: 0 == 0.0001
attach(bayescan)
class(bayescan$Q_VALUE)
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 
# Round the values
bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))
# Add a column for the type of selection grouping based on a Q-VALUE < 0.05
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.01,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 
# Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05)
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]
# Check the number of SNPs belonging to each category
xtabs(data=bayescan, ~SELECTION) 
#write.table(neutral, "neutral.txt", row.names=F, quote=F)
#write.table(balancing, "balancing.txt", row.names=F, quote=F) 
#write.table(positive, "positive.txt", row.names=F, quote=F) 

# Transformation Log of the Q value in order to plot SNPs
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 
x_title="Log(q-value)" 
y_title="Fst" 
ggplot(bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+   
  theme_classic()



###------ Make a manhattan plot, point size by allele frequency, color by adjusted pval

sites=read.table("GBRa.mafs.gz",header=T)

mh=sites[,c(1,2,5)]
mh$pos.mb=mh$position/1e+6
names(mh)[1:3]=c("chrom","pos","maf")
mh$logq=bayescan$LOG10_Q
mh$fst=bayescan$FST
mh$qval=bayescan$Q_VALUE
mh$maf3=mh$maf^3
mh=mh[1:9663,]
ggplot(mh,aes(pos.mb,logq))+
  geom_point(shape = 21, colour = "grey20", aes(size=maf3,fill=logq))+
  scale_size_continuous(breaks=c(0.2,0.4,0.6)^3,labels=c(0.2,0.4,0.6))+
  scale_fill_gradient(low="grey80",high="coral")+
  theme_bw() + 
  labs(size = "maf")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("position,Mb")

top30<-head(mh[order(mh$logq,decreasing=T),], 30)
#write.csv(top30, "top30SNPs_GBR.csv")



###------ Make a manhattan plot with chromosome distinctions
# Based on scripts from https://r-graph-gallery.com/101_Manhattan_plot.html by Yan Holtz
mh$chrom<-str_replace(mh$chrom, "chr", "")
mh$chrom<-str_replace(mh$chrom, "Sc", "")
mh$chrom<-str_replace(mh$chrom, "xf", "")
mh$chrom<-str_replace(mh$chrom, "xp", "")
mh$chrom<-as.numeric(mh$chrom)
don <- mh %>%
  group_by(chrom) %>%
  summarise(chr_len=max(pos)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  left_join(mh, ., by=c("chrom"="chrom")) %>%
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)
axisdf = don %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(don, aes(x=BPcum, y=logq)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 14)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


snps<-head(don[order(don$logp,decreasing=T),])



###------ Make the plot without numbered labels and with highlighted points
highlight_df <- don %>% 
  filter(chrom==9, pos==3997448| pos==4014530 | pos==4044872 )
signif <- don %>% 
  filter(qval<0.01)

ggplot(don, aes(x=BPcum,fst)) +
  
  # Show all points
  geom_point( aes(size=maf3, color=as.factor(chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=signif, aes(x=BPcum, y=fst), color='deeppink', size=1.3)+
  geom_point(data=highlight_df, aes(x=BPcum, y=fst), color='deeppink1', size=2.5)+
  
  
  labs(x ="Chromosome")+
  labs(y ="Fst")+
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )






