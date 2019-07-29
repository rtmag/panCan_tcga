suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

normal_sample_number <- read.table(
         pipe('grep -c "NORMAL" /root/TCGA/Rnbeads/*/*_info_withNormal.csv|perl -pe "s/\\//\t/g"|perl -pe "s/\\:/\t/g"|cut -f5,7'),
         stringsAsFactors=FALSE)

normal_sample_number <- normal_sample_number[normal_sample_number[,2]>9,1]

for(i in 1:length(normal_sample_number)){
     tcga.dir <- list.dirs(path = "/root/TCGA/Rnbeads/", full.names = TRUE, recursive = FALSE)
     command <- paste0("/root/TCGA/Rnbeads/",normal_sample_number[i],"/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")
     rnb.set.norm=load.rnb.set(command)

}

