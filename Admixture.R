setwd("~/Documents/Ds9") 
library(RColorBrewer)
# Get ID and pop info for each individual from Matz et al. 2018
site<-read.table("2bRAD_indpop.txt")
names(site)<-c("INDIVIDUALS", "SITE")
bams<-read.table("bams.nr")
names(bams)<-"INDIVIDUALS"
bams$INDIVIDUALS<- gsub(".bam","",bams$INDIVIDUALS)
pop<-merge(bams, site, by="INDIVIDUALS")


#For K=2:
# Import admixture proportions from NGSadmix:
q2<-read.table("GBR2.qopt")
# Make a barplot (ordered by population)
ord<-order(pop$INDIVIDUALS)
CSgt = colorRampPalette(colors = c('lightsteelblue2','lightsteelblue4')) # create a colorscale for plotting the genotype frequency
barplot(t(q2)[,ord],
        col=CSgt(2),
        names=site$INDIVIDUALS[ord],
        las=2,
        space=0,
        border=NA,
        ylab="Admixture proportions for K=2")



#For K=3:
# Import admixture proportions from NGSadmix:
q3<-read.table("GBR3.qopt")
# Make a barplot (ordered by population)
ord<-order(pop$INDIVIDUALS)
CSgt = colorRampPalette(colors = c('lightsteelblue4', 'lightsteelblue2', 'lightsteelblue3')) # create a colorscale for plotting the genotype frequency
barplot(t(q3)[,ord],
        col=CSgt(3),
        names=site$INDIVIDUALS[ord],
        las=2,
        space=0,
        border=NA,
        ylab="Admixture proportions for K=3")






