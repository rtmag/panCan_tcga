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
    dmc_table <- read.csv(path)
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
    # Labels
     path <- paste0("/root/TCGA/Rnbeads/",tcga_name,"/",tcga_name,"_TP53_mutation_info_withNormal.csv")
     TUMOR = read.csv(path,header=TRUE)
     TUMOR = as.character(TUMOR$Variant_Classification)
     TUMOR[TUMOR!="NORMAL"] = "TUMOR"
     track= TUMOR 
     track[track=="TUMOR"]=1
     track[track=="NORMAL"]=2
     track=as.numeric(track)
     colores=c("red","grey")
     clab=as.character(colores[track])
     colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))
    meth.norm.sig=meth.norm[abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]
    meth.norm.sig= meth.norm.sig[complete.cases(meth.norm.sig),]
    meth.norm.sig= meth.norm.sig[apply(meth.norm.sig,1,sd)>0,]
    
    path <- paste0(tcga.dir[i],"/",tcga_name,"_Tumor_VS_Normal_heatmap.png")
    png(path,width= 3.25, height= 3.25,units="in", res=1200,pointsize=4)
    heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
              labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
    dev.off()
    
    #do Volcano
    path <- paste0(tcga.dir[i],"/",tcga_name,"_Tumor_VS_Normal_volcanoPlot.pdf")
    pdf(path)
              smoothScatter(dmc_table$mean.diff,-log10(dmc_table$diffmeth.p.adj.fdr),
              ylab=expression('-Log'[10]*' Q-values'), xlab="Differential Beta-score" )
              abline(v=-.25,lty = 2,col="grey")
              abline(v=.25,lty = 2,col="grey")
              abline(h=-log10(0.05),lty = 2,col="grey")
              legend("topright", paste("NORMAL",length(which(dmc_table$mean.diff>.25 & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
              legend("topleft", paste("TUMOR",length(which(dmc_table$mean.diff<(-.25) & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
    dev.off()
}
