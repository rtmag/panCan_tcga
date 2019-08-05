suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

num.cores <- 20
parallel.setup(num.cores)

tcga.dir <- list.dirs(path = "/root/TCGA/TF_METH/", full.names = TRUE, recursive = FALSE)

for(i in 1:length(tcga.dir)){
    # get name of TCGA
    tcga_name <- gsub(".+\\/\\/","",tcga.dir[i],perl=TRUE)
    # hypo hyper cpg
    path <- paste0(tcga.dir[i],"/",tcga_name,"_dmc_table.csv")
    meth_table <- read.csv(path)
    # read rdata
    
    #do heatmap
    
    
}
