library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
setwd("~/Documents/Ds9")
load("traits.Rdata")

gain_mean <- traits %>% 
  group_by(ori, tra) %>% 
  summarize(average = mean(GAIN, na.rm=TRUE)) %>%
  ungroup()
Kgain_mean<-filter(gain_mean, ori == "K")
Ogain_mean<-filter(gain_mean, ori == "O")

ggplot(traits, aes(fct_rev(tra), GAIN, fill=fct_rev(ori)))+
  geom_boxplot()+
  scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  labs(fill='Origin site')+
  scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  ggtitle("by origin") +
  ylab("Gains")+
  xlab("Transplant site")+
  theme_bw()


f<-ggplot(Orpheus, aes(fct_rev(tra),GAIN))+
  geom_boxplot(fill="#F8766D")+
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  #geom_point(data = filter(allr.dat, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #labs(fill='genotype')+
  #scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  facet_wrap(~fct_rev(ori)) +
  #ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Transplant site")+
  theme_bw()

g<-ggplot(Keppel, aes(fct_rev(tra),GAIN))+
  geom_boxplot(fill="#00BFC4")+
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  #geom_point(data = filter(allr.dat, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #labs(fill='genotype')+
  #scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  facet_wrap(~fct_rev(ori)) +
  #ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Transplant site")+
  theme_bw()
ggarrange(f,g)




######---------------- Contrasting by genotype (using significant SNP from 2brad genotypes)
conds2<-read.csv("DS9_genotypes.csv")[,-1] # This is conds2 containing genotypes (from Gene_expression.R)
names(traits)[1]<-"sam"
gt_traits<-merge(traits, conds2, by=c("sam", "ori", "tra"))
gt_traits[gt_traits$ori == "O", "ori"] <- "Orpheus origin"
gt_traits[gt_traits$ori == "K", "ori"] <- "Keppel origin"
gt_traits[gt_traits$tra == "K", "tra"] <- "Keppel"
gt_traits[gt_traits$tra == "O", "tra"] <- "Orpheus"

ggplot(gt_traits, aes(fct_rev(genotype), GAIN, fill=fct_rev(tra)))+
  geom_boxplot()+
  geom_point(data = filter(gt_traits, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  labs(fill='Transplant site')+
  scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  facet_wrap(~fct_rev(ori)) +
  ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Genotype")+
  theme_bw()

ggplot(gt_traits, aes(fct_rev(tra), GAIN, fill=fct_rev(genotype)))+
  geom_boxplot()+
  scale_fill_manual(values=CSgt(3))+
  #geom_point(data = filter(gt_traits, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  labs(fill='genotype')+
  #scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  #facet_wrap(~fct_rev(ori)) +
  ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Transplant site")+
  theme_bw()

gt_traits[gt_traits$tra == "Orpheus", "tra"] <- "Orpheus transplant"
gt_traits[gt_traits$tra == "Keppel", "tra"] <- "Keppel transplant"

CSgt = colorRampPalette(colors = c('skyblue1', 'slateblue1', 'slateblue4')) # create a colorscale for plotting the genotype frequency
Orpheus<-gt_traits %>% dplyr::filter(ori=="Orpheus origin")
o<-ggplot(Orpheus, aes(fct_rev(genotype),GAIN, fill=fct_rev(genotype)))+
  geom_boxplot()+
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  scale_fill_manual(values=CSgt(3))+
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  #geom_point(data = filter(allr.dat, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #labs(fill='genotype')+
  #scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  facet_wrap(~fct_rev(tra)) +
  ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Orpheus origin")+
  theme_bw()

CSgt2 = colorRampPalette(colors = c('slateblue1', 'slateblue4')) # create a colorscale for plotting the genotype frequency
Keppel<-gt_traits %>% dplyr::filter(ori=="Keppel origin")
k<-ggplot(Keppel, aes(fct_rev(genotype), GAIN,fill=fct_rev(genotype)))+
  geom_boxplot()+
  #stat_summary(fun.data = give.n, geom = "text", fun = median) +
  scale_fill_manual(values=CSgt2(2))+
  #scale_x_discrete(labels=c('Orpheus', 'Keppel'))+
  #geom_point(data = filter(allr.dat, ind == "K9" | ind =="O1"), aes(fct_rev(genotype), color = fct_rev(tra)), show.legend = FALSE) +
  #labs(fill='genotype')+
  #scale_fill_discrete(labels=c('Orpheus', 'Keppel'))+
  facet_wrap(~fct_rev(tra)) +
  #ggtitle("by genotype") +
  ylab("Gains")+
  xlab("Keppel origin")+
  theme_bw()

ggarrange(o,k)
