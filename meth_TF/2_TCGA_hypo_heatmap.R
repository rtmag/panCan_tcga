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
    
    # read rdata
    path <- paste0(tcga.dir[i],"/","rnb.set.norm.filtered.RData.zip")
    rnb.set.norm=load.rnb.set(path)
    meth.norm = meth(rnb.set.norm,row.names=T)

    # hypo hyper cpg
    path <- paste0(tcga.dir[i],"/",tcga_name,"_dmc_table.csv")
    meth_table <- read.csv(path)
        #Read Annotation
        library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
        data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
        annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    
    hi_normal = meth.norm[dmc_table$mean.diff>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]
    hi_normal= hi_normal[complete.cases(hi_normal),]
    hi_normal= hi_normal[apply(hi_normal,1,sd)>0,]

    hi_tumor = meth.norm[dmc_table$mean.diff<(-.25) & dmc_table$diffmeth.p.adj.fdr<0.05,]
    hi_tumor= hi_tumor[complete.cases(hi_tumor),]
    hi_tumor= hi_tumor[apply(hi_tumor,1,sd)>0,]

    hi_tumor_anno = annotation.table[rownames(annotation.table) %in% rownames(hi_tumor),c(1,2,2)]
    path <- paste0(tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.bed")
    write.table(hi_tumor_anno, path,sep="\t",quote=F,row.names=F,col.names=F)
    hi_tumor_anno[,2] = hi_tumor_anno[,2]-50
    hi_tumor_anno[,3] = hi_tumor_anno[,3]+50
    path <- paste0(tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg_100bp.bed")
    write.table(hi_tumor_anno, path,sep="\t",quote=F,row.names=F,col.names=F)

    hi_normal_anno = annotation.table[rownames(annotation.table) %in% rownames(hi_normal),c(1,2,2)]
    path <- paste0(tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.bed")
    write.table(hi_normal_anno, path,sep="\t",quote=F,row.names=F,col.names=F)
    hi_normal_anno[,2] = hi_normal_anno[,2]-50
    hi_normal_anno[,3] = hi_normal_anno[,3]+50
    path <- paste0(tcga.dir[i],"/",tcga_name,"_hi_normal_cpg_100bp.bed")
    write.table(hi_normal_anno, path,sep="\t",quote=F,row.names=F,col.names=F)

    
    
    #do heatmap
    
    
}
