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
################################################################################################
################################################################################################
################################################################################################
tcga.dir <- list.dirs(path = "/root/TCGA/TF_METH/", full.names = TRUE, recursive = FALSE)

for(i in 1:length(tcga.dir)){
    tcga_name <- gsub(".+\\/\\/","",tcga.dir[i],perl=TRUE)

command <- paste0(tcga.dir[i],"/",tcga_name,"_DMcpg_annStats_pie.pdf")
pdf(command)
par(mfrow=c(1,2))
command <- paste0("cut -f1,2,4 ",tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.annStats")
res=read.table(pipe(command), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdownfin = tdown[round(tdown/sum(tdown)*100,digits=2)>=2]
other = sum(tdown[round(tdown/sum(tdown)*100,digits=2)<2])
names(other) = paste0("Other ",round(other/sum(tdown)*100,digits=2),"%")
tdownfin = c(tdownfin,other)
pie(sort(tdownfin), main="Normal",cex=.8)

command <- paste0("cut -f1,2,4 ",tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.annStats")
res=read.table(pipe(command), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdownfin = tdown[round(tdown/sum(tdown)*100,digits=2)>=2]
other = sum(tdown[round(tdown/sum(tdown)*100,digits=2)<2])
names(other) = paste0("Other ",round(other/sum(tdown)*100,digits=2),"%")
tdownfin = c(tdownfin,other)
pie(sort(tdownfin), main="Tumor",cex=.8)
dev.off()
####################################
command <- paste0(tcga.dir[i],"/",tcga_name,"_DMcpg_annStats_barplot_enrichment.pdf")
pdf(command)
par(mfrow=c(1,2))

command <- paste0("cut -f1,2,4 ",tcga.dir[i],"/",tcga_name,"_hi_normal_cpg.annStats")
res1=read.table(pipe(command), sep="\t",header=F)
i1 = which(res1[,1]=="Annotation")[2]+1
i2 = dim(res1)[1]
res1 = res1[ i1:i2,]
tdown_th = as.numeric(as.character(res1[,2]))
tdownfin_th = tdown_th[round(tdown_th/sum(tdown_th)*100,digits=2)>=2]
tdown = as.numeric(as.character(res1[,3]))
names(tdown) = res1[,1]
tdown1 = tdown[round(tdown_th/sum(tdown_th)*100,digits=2)>=2]

command <- paste0("cut -f1,2,4 ",tcga.dir[i],"/",tcga_name,"_hi_tumor_cpg.annStats")
res2=read.table(pipe(command), sep="\t",header=F)
i1 = which(res2[,1]=="Annotation")[2]+1
i2 = dim(res2)[1]
res2 = res2[ i1:i2,]
tdown = as.numeric(as.character(res2[,3]))
names(tdown) = res2[,1]
tdown2 = tdown[round(tdown_th/sum(tdown_th)*100,digits=2)>=2]

barplot(sort(tdown1),las=2,ylim=c(min(c(tdown1,tdown2)),max(c(tdown1,tdown2))),ylab="Log2 Enrichment",col="lightblue3",main="Normal")
abline(h=0)
barplot(sort(tdown2),las=2,ylim=c(min(c(tdown1,tdown2)),max(c(tdown1,tdown2))),ylab="Log2 Enrichment",col="lightblue3",main="Tumor")
abline(h=0)
dev.off()    
}
################################################################################################
################################################################################################
################################################################################################
