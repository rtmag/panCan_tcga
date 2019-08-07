options(bitmapType="cairo")
options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

tcga.dir <- list.dirs(path = "/root/TCGA/TF_METH/", full.names = TRUE, recursive = FALSE)


for(i in 1:length(tcga.dir)){
    # get name of TCGA
    tcga_name <- gsub(".+\\/\\/","",tcga.dir[i],perl=TRUE)
  
    # annotate hi norm cpg
    command <- paste0("awk '{print $1\"\t\"$2-1\"\t\"$3}' ",
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.bed|"),
                      "annotatePeaks.pl - hg19 -annStats ",
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.annStats > "),
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.anno")
                     )        
    system(command)

    # annotate hi tumor cpg
    command <- paste0("awk '{print $1\"\t\"$2-1\"\t\"$3}' ",
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.bed|"),
                      "annotatePeaks.pl - hg19 -annStats ",
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.annStats > "),
                      paste0(tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.anno")
                     )        
    system(command)
  }
