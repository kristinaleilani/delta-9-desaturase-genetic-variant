# All sequence processing was run on Lonestar6 at the Texas Advanced Computing Center

# --------- Installations ----------
# download 2bRAD processing scripts:
git clone https://github.com/z0on/2bRAD_denovo.git

# download tools to download SRA datasets from NCBI:
# esearch : see here for installation (on TACC):
https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
cd
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
# SRA toolkit (on TACC, pick the one for ubuntu 64):
cd
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar vxf sratoolkit.current-ubuntu64.tar.gz
cd
nano .bashrc
# section 2:
   export PATH=$HOME/sratoolkit.3.0.1-ubuntu64/bin:$PATH
   export PATH=$HOME/edirect:$PATH
# save (Ctl-O, Ctl-X)

# All other packages installed via conda





#---------------------------------- 2bRAD data analysis --------
# Based on scripts from: https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh

# Download fastq from Matz et al. 2018 paper:
export BioProject=PRJNA434194
$HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv

>gets
for A in `cat $BioProject.SRR`;do 
echo "fastq-dump-orig.3.0.0 $A">>gets;
done
ls6_launcher_creator.py -j gets -n gets -a IBN21018 -e kblack@utexas.edu -t 20:00:00 -w 12 -q normal
sbatch gets.slurm

#import A. millepora reference genome (from Fuller et al. 2020)
cdw
cd db
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna.gz 

# Also, import symbiont reference genomes (stored in a shared directory on TACC) and concatenate to coral genome
cd /work/06909/cbscott/sym_references/
cp Amil_symABCD* /scratch/07090/kblack/Ds9/split/OK_concatenated
export GENOME_FASTA=Amil_symABCD.fasta
samtools faidx $GENOME_FASTA
grep ">chr" Amil.v2.01.chrs.fasta # Chromosomes 20-24 belong to symbionts A-D, respectively

# trimming and subsampling to 3M filtered reads max
export TagLen=36  
export MatchFrac=0.95 # up to 5% divergence 
>trim2
for file in *.fastq; do
echo "cutadapt -q 15,15 -a AGATCGGA  -m $TagLen -l $TagLen -o ${file/.fastq/}.trim0 $file > ${file}_trimlog.txt && head -12000000 ${file/.fastq/}.trim0 > ${file/.fastq/}.trim && rm ${file/.fastq/}.trim0" >> trim2; done
ls6_launcher_creator.py -j trim2 -n trim2 -a IBN21018 -e kblack@utexas.edu -t 1:00:00 -w 48 -q normal
sbatch trim2.slurm

# mapping, converting to bams, indexing
export GENOME_FASTA=Amil_symABCD.fasta
>maps
for F in `ls *.trim`; do
echo "bowtie2 --local --no-unal -x $GENOME_FASTA -U ${F} -S ${F/.trim/}.sam && samtools sort -O bam -o ${F/.trim/}.bam ${F/.trim/}.sam && samtools index ${F/.trim/}.bam">>maps
done
ls6_launcher_creator.py -j maps -n maps -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 12 -q normal
sbatch maps.slurm

mkdir split
cp *bam* split

# Split out coral reads to keep (remove reads mapping to symbiont genomes)
cd /work/06909/cbscott/sym_references/
cp Amil.v2.01.chrs.fasta /scratch/07090/kblack/Ds9/split/OK_concatenated
cd /scratch/07090/kblack/Ds9/split/OK_concatenated
export GENOME=Amil.v2.01.chrs.fasta
samtools faidx $GENOME
export GENOME_FAI=Amil.v2.01.chrs.fasta.fai
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $GENOME_FAI > Amil.bed

cp Amil.bed split

cd split
>split
for F in `ls *.bam`; do
echo "samtools view -L Amil.bed -o ${F/.bam/}_Amil.bam $F" >>split
done
ls6_launcher_creator.py -j split -n split -a IBN21018 -e kblack@utexas.edu -t 00:30:00 
sbatch split.slurm

# Rename samples
awk -F'\t' 'system("mv " $1 " " $2)' SampleID.txt
# Reindex new bams names
>maps2
for F in `ls *.bam`; do
echo "samtools index $F">>maps2
done
ls6_launcher_creator.py -j maps2 -n maps2 -a IBN21018 -e kblack@utexas.edu -t 0:10:00  
sbatch maps2.slurm

# Or remove _Amil from all file names
mkdir rename
cp *_Amil.bam rename
cd rename
for filename in *.bam; do 
    [ -f "$filename" ] || continue
    mv "$filename" "${filename//_Amil/}"

done



# quality assessment, removing bams with log(coverage)<3SD
# also imposing genotyping rate cutoff - 50%
conda activate base # Use angsd installed locally (conda installed angsd doesn't work sometimes)
module load Rstats
export MinIndPerc=0.8
FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 30"
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo "ls *.bam > bams && angsd -b bams -r chr1 -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R dd >qualRanks">a0
ls6_launcher_creator.py -j a0 -n a0 -a IBN21018 -e kblack@utexas.edu -t 1:00:00  
sbatch a0.slurm
# look at quality of reads in dd.pdf

# Detecting and removing clones (see hctree.pdf and resulting bams.nr)
export MinIndPerc=0.8
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.qc | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.qc -GL 1 $FILTERS0 $TODO0 -P 12 -out GBR && Rscript ~/bin/detect_clones.R bams.qc GBR.ibsMat 0.15">a1
ls6_launcher_creator.py -j a1 -n a1 -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 1
sbatch a1.slurm
# retained 118 out of 126 bams

# Prepping beagle file for admixture (low quality bams and clones removed)
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'cat bams.nr | sort > bams.NR && mv bams.NR bams.nr && export NIND2=`cat bams.nr | wc -l`; export MI2=`echo "($NIND2*$MinIndPerc+0.5)/1" | bc`' >calc2
echo "source calc2 && angsd -b bams.nr -GL 1 $FILTERS1 $TODO1 -P 12 -out GBR2 && Rscript ~/bin/pcaStructure.R GBR2.ibsMat > pcaStruc.txt">a2
ls6_launcher_creator.py -j a2 -n a2 -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 1 
sbatch a2.slurm




#### --------- Admixture ------------------------
# ngsAdmix makes 3 qopt files for 1-3 pops
echo 'for K in `seq 3` ; do  NGSadmix -likes GBR2.beagle.gz -K $K -P 12 -o GBR${K}; done' >adm3
ls6_launcher_creator.py -j adm3 -n adm3 -a IBN21018 -e kblack@utexas.edu -t 1:00:00 -w 1 
sbatch adm3.slurm





#### -------- Make  VCF for Sambada ------------------
echo "bcftools mpileup -Ou -f Amil.v2.01.chrs.fasta -b bams.nr | bcftools call -mv -Ob -o GBRall.bcf">vcf2
ls6_launcher_creator.py -j vcf2 -n vcf2 -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 24 
sbatch vcf2.slurm
# 569263 sites
bcftools view GBRall.bcf > GBRall.vcf #Can run the rest on ls6 on idev
bcftools view -i 'F_MISSING < 0.1' GBRall.vcf -Ov -o GBR1all.vcf # Remove variants which have missing values less than 10%.
# 56047 sites
bcftools view -i 'MAF > 0.05' GBR1all.vcf -Ov -o GBR2all.vcf # Remove by minor allele frequency.
# 21687 sites
bcftools view -c 1 GBR2all.vcf -o GBRall.vcf -Oz # Remove monomorphic sites
# 21688 sites

# Set minimum quality
echo "bcftools filter -i'%QUAL>50' GBRall.vcf > GBRq.vcf">filter1
ls6_launcher_creator.py -j filt -n filt -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 24 -q normal
sbatch filt.slurm
# 21558 sites

# Filter out indels
echo "vcfutils.pl varFilter GBRq.vcf > GBRqv.vcf" > filt1
ls6_launcher_creator.py -j filt1 -n filt1 -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 24 -q normal
sbatch filt1.slurm
# 20993 sites

# Set ID column names
echo "bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' GBRqv.vcf > GBRqvid.vcf">ann
ls6_launcher_creator.py -j ann -n ann -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 24 -q normal
sbatch ann.slurm





#### --------- Sambada------------------------
echo "/home1/07090/kblack/sambada-0.8.3-ubuntu/binaries/sambada GBR_OK_allSNPs_nopops-param.txt GBR_OK_allSNPs_nopops-environmental-matrix.txt GBR_OK_allSNPs-genotype-matrix.txt" > sam
ls6_launcher_creator.py -j sam -n sam -a IBN21018 -e kblack@utexas.edu -t 10:00:00 -w 24 -q normal
sbatch sam.slurm





####### --------- Bayescan: identifying Fst outliers ---------------------------------------
# Running ANGSD on all bams with filters that do not disturb AFS, to find sites to work on.
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11"
echo "angsd -b bams.nr -GL 1 -P 12 $FILTERS $TODO -minInd 100 -out oksites">ac
ls6_launcher_creator.py -j ac -n ac -t 0:30:00 -e kblack@utexas.edu -w 1 -a IBN21018 
sbatch ac.slurm

# collecting sites that are not likely to be lumped paralogs (i.e., with majority of calls being heterozygotes)
idev
zcat oksites.geno.gz | python ~/bin/HetMajorityProb.py | awk "\$6 < 0.75 {print \$1\"\t\"\$2}" > goodsites
angsd sites index goodsites
exit

# making VCF for bayescan (only variable sites):
export FILTERS1="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 100 -snp_pval 1e-5 -minMaf 0.05 "
export TODO1="-gl 1 -dopost 1 -domajorminor 1 -domaf 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -docounts 1"
echo "angsd -b bams.nr -sites goodsites $FILTERS1 $TODO1 -P 1 -out GBRa" > makebcf
ls6_launcher_creator.py -j makebcf -n makebcf -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 1
sbatch makebcf.slurm

conda activate bcftools
bcftools view GBRa.bcf > GBRa.vcf
grep -v "#" GBRa.vcf | cut -f 1,2 > bsa.sites #11005 sites 

# creating file with population designations
grep "#CHR" GBRa.vcf | cut -f 10- | perl -pe 's/\s/\n/g' >inds
cat bams.nr | perl -pe 's/^((.).+)/$2/' >pops
paste inds pops -d '\t' >bspops

# create a file called vcf2bayescan.spid :
echo "############
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan.spid

# Install PGDspider (conda doesn't work)
cd ~/bin
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.1.5.zip
unzip PGDSpider_2.1.1.5.zip
cd -
# Install Bayescan (conda doesn't work)
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip 
cd BayeScan2.1/binaries
chmod u+x BayeScan2.1_linux64bits

# converting vcf to bayescan format
echo "java -jar ~/bin/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile GBRa.vcf  -outputfile GBRa.bayescan -spid vcf2bayescan.spid " >bconv
ls6_launcher_creator.py -j bconv -n bconv -t 0:15:00 -a IBN21018 -e kblack@utexas.edu -w 1 
sbatch bconv.slurm

# Prior odds =10 (for every 10 neutral loci in the data set, odds are that 1 locus is under selection)
echo "~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits GBRa.bayescan -threads=24">bsc
ls6_launcher_creator.py -j bsc -n bsc -t 24:00:00 -a IBN21018 -e kblack@utexas.edu -w 1 -q normal
sbatch bsc.slurm

# scp sites (GBRa.mafs.gz) and bayescan results (GBRa.baye_fst.txt)
# use Bayescan_plots.R to examine results 





####### --------- Gene expression ---------------------------------------

# First download Tagseq data from Dixon 2018 paper and map to coral genome
# Need sam files:
>maps
for F in `ls *.fastq`; do
echo "samtools sort -O bam -o ${F/.fastq/}.bam ${F/.fastq/}.sam">>maps
done
ls6_launcher_creator.py -j maps -n maps -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 12 -q normal 
sbatch maps.slurm

# Download metadata
cd /scratch/07090/kblack/Ds9/rna
export BioProject=PRJNA265072
$HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo > $BioProject.fullMeta.csv
# Make list of samples (SampleID_rna.txt) and import 
cat -A SampleID_rna.txt | head # See if there are special invisible characters
LC_ALL=C sed 's/[^[:blank:][:print:]]//g' SampleID_rna.txt > SampleID_rna.new # Remove them
sed 's/\^I//g' SampleID_rna.new > SampleID_rna.new2 
# Rename samples
awk -F'\t' 'system("mv " $1 " " $2)' SampleID_rna.new2

# Reindex new bams names
>maps2
for F in `ls *.bam`; do
echo "samtools index -c $F">>maps2
done
ls6_launcher_creator.py -j maps2 -n maps2 -a IBN21018 -e kblack@utexas.edu -t 0:10:00  
sbatch maps2.slurm

### Download Amil annotations
cd $STOCKYARD 
mkdir Amil
cd Amil
wget https://www.dropbox.com/sh/h75426cew6h1ttc/AABI672CSygMXY1enEX1jwfCa/._Amil.coding.gff3?dl=1

# get gff
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.gff.gz .  
gunzip GCF_013753865.1_Amil_v2.1_genomic.gff.gz
sed -i -e 's/NC_058066.1/chr1/g' *.gff
sed -i -e 's/NC_058067.1/chr2/g' *.gff
sed -i -e 's/NC_058068.1/chr3/g' *.gff
sed -i -e 's/NC_058069.1/chr4/g' *.gff
sed -i -e 's/NC_058070.1/chr5/g' *.gff
sed -i -e 's/NC_058071.1/chr6/g' *.gff
sed -i -e 's/NC_058072.1/chr7/g' *.gff
sed -i -e 's/NC_058073.1/chr8/g' *.gff
sed -i -e 's/NC_058074.1/chr9/g' *.gff
sed -i -e 's/NC_058075.1/chr10/g' *.gff
sed -i -e 's/NC_058076.1/chr11/g' *.gff
sed -i -e 's/NC_058077.1/chr12/g' *.gff
sed -i -e 's/NC_058078.1/chr13/g' *.gff
sed -i -e 's/NC_058079.1/chr14/g' *.gff
# Remove bad gene
grep 'LOC114946994' GCF_013753865.1_Amil_v2.1_genomic.gff
sed -i '/LOC114946994/d' GCF_013753865.1_Amil_v2.1_genomic.gff

#run featurecounts #
>runFeatureCounts
for file in *.bam
do echo "featureCounts -a GCF_013753865.1_Amil_v2.1_genomic.gff -p -t gene -g gene -o feature_counts_out.tsv -T 48 --primary *.bam" >> runFeatureCounts
done
ls6_launcher_creator.py -n runFeatureCounts -j runFeatureCounts -w 48 -N 4 -a IBN21018 -e kblack@utexas.edu -t 10:00:00 -q normal
sbatch runFeatureCounts.slurm

sed -i 1d feature_counts_out.tsv

# scp feature_counts_out.tsv to explore in DESeq2









#####---------- Haplotype selection stats ---------
# based on scripts from: https://garud.eeb.ucla.edu/selection-scans/

# import chr9_onepop.txt and chr9.vcf (from Fuller et al. 2020)
lassip --vcf chr9.vcf --avg-spec --hapstats --salti --winsize 1000 --winstep 25 --out chr9_salti --pop chr9_onepop.txt
# Change column names
sed -e '1s/GBR_h12/H12/' -e '1s/GBR_h2h1/H2H1/' chr9.GBR.lassip.hap.stats > chr9_H12

# Using SelectionHapStats
#Calculate H12
python2 ~/SelectionHapStats/scripts/H12_H2H1.py chr9_H12_input -w 1000 88 -o chr9_H12_output
# Find peaks
python2 ~/SelectionHapStats/scripts/H12peakFinder.py chr9_H12_output -o chr9_H12_peaks -t 0.02
# Manhattan plot
Rscript ~/SelectionHapStats/scripts/H12_viz.R chr9_H12_output chr9_H12_peaks chr9_H121000scan.pdf 10
# SFS
Rscript ~/SelectionHapStats/scripts/hapSpectrum_viz.R chr9_H12_peaks hapSpectrum.pdf 10 88
# Visualize chromosome, Focus on center of significant peak near Ds9
bash ~/SelectionHapStats/scripts/visualizeGenomicData.sh chr9_H12_input chr9_H12_peaks 4063469 401 88 chr9_genomicViz4063469.pdf







