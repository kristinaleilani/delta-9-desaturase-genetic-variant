setwd("~/Documents/Ds9")
#install.packages("chromoMap")
library(chromoMap)

# chromosome files
chr_file_1 = "chromosome_file.txt"
head(read.table(chr_file_1,sep = "\t"))

chr_file_2 = "chromzoom_file.txt"
head(read.table(chr_file_2,sep = "\t"))


# annotation files
anno_file_1 = "annotation_file.txt"
head(read.table(anno_file_1,sep = "\t"))
anno_file_1 = "annotations2.txt"
head(read.table(anno_file_1,sep = "\t"))


library(chromoMap)
chromoMap(chr_file_1,anno_file_1)

chromoMap("chromosome_file.txt","annotation_file.txt",
          labels=T,
          chr_color = c("lightgrey"),
          anno_col = c("deeppink"),
          label_angle = -65,
          chr_length = 6,
          chr_width = 25,
          canvas_width = 800)

chromoMap("chromzoom_file.txt","annotations2.txt",
          labels=T,
          chr_color = c("lightgrey"),
          anno_col = c("deeppink"),
          label_angle = -65,
          chr_length = 6,
          chr_width = 25,
          canvas_width = 800)

chromoMap("chromzoom_file.txt","annotation_file.txt",
          labels=T,
          chr_color = c("lightgrey"),
          anno_col = c("deeppink"),
          #data_based_color_map = T,
          #data_type = "numeric",
          label_angle = -65,
          chr_length = 6,
          chr_width = 25,
          canvas_width = 800,
          vertical_grid = T,
          grid_array = c(3997448, 4014530, 4034118, 4034119, 4044872, 4063469, 4068048, 4073718, 4073719, 4073720),
          grid_text = c("SNP1","SNP2","SNP3","SNP4","SNP5", "SNP6", "H12", "SNP8", "SNP9", "SNP10"))
