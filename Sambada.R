### This script is based on the Seascape Genomics Course (Physalia 2021), and this section was taught by Oliver Selmoni

### Prepare input files for SamBada
setwd("~/Documents/Ds9") 
library('vcfR') 
library('poppr') 
library('visreg') 
library('qvalue') 
library("raster")
library("rnaturalearth") 
library("rnaturalearthdata") 
library("ggplot2")
library("ggsn")
library("corrplot")
library('scatterpie')
library('dplyr')
library('maptools')
source('Sambada_functions.R')

###-------------------------------- Prepare the genetic input file
#Load the vcf 
VCF=read.vcfR('GBRall.vcf') 
# Transform to a genind object
genind <- vcfR2genind(VCF) # 119 individuals; 19,691 loci
# Remove loci with too many missing
genind.loci <- missingno(genind, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE) # No loci with missing values above 10% found.
# Remove individuals with too many missing
genind.loci.ind <- missingno(genind.loci, type = "geno", cutoff = 0.05, quiet = FALSE, freq = FALSE) # 7 genotypes contained missing values greater than 5%
# Apply a MAF criteria
genind.ind.loci.maf <- informloci(genind.loci.ind, cutoff = 2/nInd(genind.loci.ind), MAF = 0.05, quiet = FALSE) # 112 individuals; 19,359 loci


# Extract the genotype matrix from the the filtered GenInd object: 
alleleTAB=genind.ind.loci.maf$tab # extract the table of allele frequencies
#alleleTAB=genind$tab 
SNPS = alleleTAB[,rep(c(T,F), times=ncol(alleleTAB)/2)] # keep only one column per SNPs (one allele frequency)
colnames(SNPS) = substr(colnames(SNPS), 1 , nchar(colnames(SNPS))-2) # reformat colnames to contain SNPid
# Check the SNP table: 
SNPS[1:10,1:10] # each SNP is encoded in 0, 1, 2 format.
# Then convert:
GT_0 <- (SNPS==0)+0 # transform GT==0 to 1, the rest to 0
GT_1 <- (SNPS==1)+0 # transform GT==1 to 1, the rest to 0
GT_2 <- (SNPS==2)+0 # transform GT==2 to 1, the rest to 0
colnames(GT_0) <- paste0(colnames(GT_0),'_0') # add GT identifier to each SNP name
colnames(GT_1) <- paste0(colnames(GT_1),'_1') # add GT identifier to each SNP name
colnames(GT_2) <- paste0(colnames(GT_2),'_2') # add GT identifier to each SNP name

# Create the genotype matrix for SamBada, with the identifiers of individuals in the first column and the remaining columns encoding the presence (1) or absence (0) of each genotype. 
SAMGT = cbind('SampleID'=rownames(SNPS), GT_0, GT_1, GT_2) 
dim(SAMGT) # 112 59683
dim(SNPS) # 112 19894

# Write the genetic input for SamBada in a file. 
write.table(SAMGT, 'GBR_OK_allSNPs-genotype-matrix.txt', sep=' ', row.names = F, col.names = T, quote=F)


###-------------------------------- Prepare the environmental input file
# First, load the environmental matrix of sampling sites.
env=read.table("GBRall_indpop.txt")
names(env)<- c("Sample", "Site")
#Orpheus coords= -18.621843537170466, 146.48303743432103
#Keppels coords= -23.170998266481792, 150.94531935664864
#Magnetic coords= -19.157070136040897, 146.79688039654795
#Wilkie coords= -13.773463561154893, 143.63885892176847
#Sudbury coords= -16.97369827771356, 146.1633153941376
env$Latitude <- ifelse(env$Site=="O", -18.6218,
                       ifelse(env$Site=="K", -23.1709,
                              ifelse(env$Site=="M", -19.1570,
                                     ifelse(env$Site=="W", -13.7734,
                                            -16.9736))))
env$Longitude <- ifelse(env$Site=="O", 146.4830,
                        ifelse(env$Site=="K", 150.9453,
                               ifelse(env$Site=="M", 146.7968,
                                      ifelse(env$Site=="W", 143.6388,
                                             146.1633))))

goodbams<-data.frame(rownames(SNPS))
names(goodbams)[1]<- "Sample"
goodbams$Sample <- gsub(".bam","",goodbams$Sample)

# Load population structure
# Import admixture
admix0=read.csv("GBR_admix.csv")[,2:5]
names(admix0)[4]<-"Sample"
admix1 <- merge(goodbams, admix0, by="Sample")

# Combine environmental table
ENV <- left_join(goodbams, env, by="Sample")
rownames(ENV) <- ENV[,1]
names(ENV)[1]<- "SampleID"
dim(ENV) #112   4
ENV$Site<-NULL
ENV$SampleID = rownames(SNPS)
ENV$pop1 = admix1$V1
ENV$pop2 = admix1$V2
ENV$pop3 = admix1$V3

# Write the environmental input for SamBada in a file. 
write.table(ENV, 'GBR_OK_allSNPs_nopops-environmental-matrix.txt', sep=' ', row.names = F, col.names = T, quote=F)



###-------------------------------- Prepare the parameter input file
paramtext = c('HEADERS YES', # ...indicates that the genotype and environmental matrix have headers. 
              'WORDDELIM " "', # ...indicates that columns are separated by space character (" ")
              paste0('NUMINDIV ', nrow(ENV)), # ... indicates the number of individuals (must be identical and in the same order between environmental and genotype matrix)
              paste0('NUMVARENV ', ncol(ENV)), # ... indicates the number of columns in the environmental matrix
              paste0('NUMMARK ', ncol(SAMGT)), # ... indicates the number of columns in the genotype matrix
              'IDINDIV SampleID', # ... indicates that the columns called "SampleID" contains the samples ids. 
              paste0('DIMMAX ', sum(substr(colnames(ENV), 1, 3) =='pop')+1 ), # ... indicates the number of dimensions for the models that we want to compute. This will be further explained later during this exercise, but the idea is to set the number of dimensions equal to the number of variables describing the population structure (here 2) + 1. 
              'POPULATIONVAR LAST', # ... indicates that the variables describing the population structure are in the last columns of the environmental table.
              'SAVETYPE END ALL' # ... indicates that SamBada should save all the models, and sort them by significance at the end of the computations. 
)

# Save the parameter file: 
write.table(paramtext, 'GBR_OK_allSNPs_nopops-param.txt', col.names=F, row.names=F, quote=F)



###--------------------------------  Run Sambada using all three files (locally or on cluster). Line to do this is pasted below:
# sambada GBR_OK_allSNPs_nopops-param.txt GBR_OK_allSNPs_nopops-environmental-matrix.txt GBR_OK_allSNPs-genotype-matrix.txt








###-------------------------------- Post processing
# load output from SamBada
samout = read.table('GBR_OK_allSNPs-genotype-matrix-Out-1.txt', header=T, stringsAsFactors = F) 

# Compute p-values and correct for false discoveries (Gscore and Waldscore)
samout$pvalueG <- 1-pchisq(samout$GscorePop, df=1)
samout$pvalueW <- 1-pchisq(samout$WaldScorePop, df=1)

# look at the distribution of p-values
hist(samout$pvalueG, breaks=1000, main='P-values G-test' )
hist(samout$pvalueW, breaks=1000, main='P-values Wald-test' )

# p-value correction:
# create the two columns that will contain the q-values
samout$qG = rep(NA, length=nrow(samout)) # qvalue Gscore
samout$qW = rep(NA, length=nrow(samout)) # qvalue Wald
for (e in unique(samout$Env_1)) { 
  
  samout$qG[samout$Env_1==e] = qvalue(samout$pvalueG[samout$Env_1==e])$qvalue # calculate qvalue G-score
  samout$qW[samout$Env_1==e] = qvalue(samout$pvalueW[samout$Env_1==e])$qvalue # calculate qvalue Wald
  
}

# filter significant SNPs at a cut-off of q<0.05.
SIGN = samout[which(samout$qG<0.05),c(1:4,17:21)]
dim(SIGN) # 252 significant models
SIGN # many refer to different genotypes of the same SNP

# focus only on the strongest association for every SNP.  
snpids = substr(SIGN$Marker, 1, nchar(SIGN$Marker)-2)
# for every group, only retain the association model with the lowest q-value of G-score (usually p-values of G-score and Wald-score are highly correlated)
sSIGN_list = (by(SIGN, snpids, function(x) { x[which.min(x$qG),]  }))
sSIGN = data.frame(t(sapply(sSIGN_list, c)))
# order according to q-value
sSIGN = sSIGN[order(unlist(sSIGN$qG)),]
sSIGN <- as.data.frame(apply(sSIGN,2,as.character))
dim(sSIGN) # 106 unique SNPs among the significant associations. 

Lon_models <- filter(sSIGN, Env_1 == "Longitude") # 32 models
Lat_models <- filter(sSIGN, Env_1 == "Latitude") # 74 models





###-------------------------------- Visualize a significant genotype-environment association
CSgt = colorRampPalette(colors = c('slateblue4', 'slateblue1', 'skyblue1')) 
gt = SAMGT[,'chr9_3997448_T_C_2'] # delta-9 desaturase SNP
gt=as.numeric(gt)
env = ENV[,'Longitude']
# Build a logistic model with glm() function setting the "binomial" family option.
MOD = glm(gt~env, family = 'binomial')
# visualize the regression
visreg(MOD, scale='response', xlab='Latitude', ylab='GT frequency')


# visualize the genotype association with the environmental variable at the SNP level.
GT=extract.gt(VCF) # extract the genotypes (gt) from the vcf object. 
ENV$SampleID = paste(ENV$SampleID,"bam",sep = ".")
SNP_gts = GT['chr9_3997448_T_C',ENV$SampleID] # gts associated with the significant SNP of interest
SNP_gts = unlist(as.vector(SNP_gts))
# plot the distribution of genotype across the values of Environmental variation
par(mfrow=c(1,2))
boxplot(ENV$Latitude~SNP_gts, xlab='Genotype of SNP: chr9_3997448_T_C_2', ylab='Longitude', col=CSgt(3))
boxplot(ENV$Longitude~SNP_gts, xlab='Genotype of SNP: chr9_3997448_T_C_2', ylab='Longitude', col=CSgt(3))


# See which base pairs are at sites within delta-9 desaturase gene region
dp2<-extract.gt(VCF, element = "GT", return.alleles = TRUE)
d2<-as.data.frame(dp2)
d2$locus<-row.names(d2)
d2[1:4,1:3]
d2[d2$locus == "chr9_3997448_T_C", ] # Keppel is mostly C
d2[d2$locus == "chr9_4014530_T_G", ] # Keppel is all G
d2[d2$locus == "chr9_4034118_G_T", ] # Keppel is all G
d2[d2$locus == "chr9_4034119_T_C", ] # Keppel is all T


###--------------------------------  Visualize the spatial distribution of a genotype-environment association. 
# Load the shoreline (download from www.ngdc.noaa.gov/mgg/shorelines)
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(142, 152), ylim = c(-24, -11)) %>%
  fortify()

env=read.table("GBRall_indpop.txt")
names(env)<- c("site_id", "Pop")
env$site_id <- gsub("a","",env$site_id)

SAMGT<-as.data.frame(SAMGT)
site_id = SAMGT$SampleID  

coords <- ENV[,3:4]

SNP_gts1<-as.data.frame(cbind(site_id,SNP_gts))
SNP_gts1 <- SNP_gts1[!is.na(SNP_gts1)]
coords1<-as.data.frame(cbind(site_id,coords,SNP_gts))
coords1 <- na.omit(coords1)

# Formatting a table to plot piecharts on the map of Australia
d<-left_join(coords1, env, by="site_id")
d0<-d[,-1]
d2<-aggregate(list(numdup=rep(1,nrow(d0))), d0, length)
d3<-d2[,3:5]
d4<-reshape(d3, idvar = "Pop", timevar = "SNP_gts", direction = "wide")
d4[is.na(d4)] <- 0
d5<-d0
d5$SNP_gts<-NULL
d5$site_id<-NULL
d6<-aggregate(list(numdup=rep(1,nrow(d5))), d5, length)
d7<-d6[,1:3]
d8<-data.frame(left_join(d4, d7, by="Pop"))
names(d8)[2:4]<-c("CC", "CT", "TT")
d8$radius <- d8$CC + d8$CT + d8$TT
d8$radius <- sqrt(sqrt(sqrt(sqrt(sqrt(d8$radius)))))
CSgt = colorRampPalette(colors = c('slateblue4', 'slateblue1', 'skyblue1')) # create a colorscale for plotting the genotype frequency
ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data = coords1, aes(x = Longitude, y = Latitude, colour=SNP_gts, size="0.5"), position=position_jitter(h=0.3, w=0.3)) + 
  guides(size="none") + 
  scale_colour_manual(values=CSgt(3), name='GT frequency')+
  theme(panel.background = element_rect(fill = "white", colour = "black"),  axis.title = element_blank()) 


p <- ggplot() +
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+  
  theme(panel.background = element_rect(fill = "white", colour = "black"),  axis.title = element_blank()) 
p + geom_scatterpie(aes(x=Longitude, y=Latitude, 
                        group=Pop, r=radius/2),
                    data=d8, cols=colnames(d8[2:4]), color=NA, alpha=.8) + 
  scale_fill_manual(values=c('slateblue4', 'slateblue1', 'skyblue1')) +
  labs(fill='Genotype')+
  coord_fixed()




