# Filters, differential CPG between tumor and normal
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
library(graphics)

rnb.set.norm=load.rnb.set("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")

TUMOR = read.csv("/root/TCGA/Rnbeads/COAD/COAD_TP53_mutation_info_withNormal.csv",header=TRUE)
TUMOR = as.character(TUMOR$Variant_Classification)
TUMOR[TUMOR!="NORMAL"] = "TUMOR"
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Tumor = TUMOR)

num.cores <- 20
parallel.setup(num.cores)
dmc <- rnb.execute.computeDiffMeth(rnb.set.norm,pheno.cols=c("Tumor"))
comparison <- get.comparisons(dmc)[1]
dmc_table <-get.table(dmc, comparison, "sites", return.data.frame=TRUE)
#
table(abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05)
library(gplots)
library(factoextra)
library(RColorBrewer)

meth.norm = meth(rnb.set.norm,row.names=T)
meth.norm.sig=meth.norm[abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]
dim(meth.norm.sig)
meth.norm.sig= meth.norm.sig[complete.cases(meth.norm.sig),]
dim(meth.norm.sig)
meth.norm.sig= meth.norm.sig[apply(meth.norm.sig,1,sd)>0,]
dim(meth.norm.sig)

track= TUMOR 
track[track=="TUMOR"]=1
track[track=="NORMAL"]=2
track=as.numeric(track)
colores=c("red","white")
clab=as.character(colores[track])

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

png("heatmap_diff_TUMOR_VS_NORMAL_FDR5p.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TUMOR","NORMAL"),fill=c("red","white"), border=T, bty="n" )
dev.off()

pdf("volcano_diff_TUMOR_VS_NORMAL.pdf")
smoothScatter(dmc_table$mean.diff,-log10(dmc_table$diffmeth.p.adj.fdr),
              ylab=expression('-Log'[10]*' Q-values'), xlab="Differential Beta-score" )
  
abline(v=-.25,lty = 2,col="grey")
abline(v=.25,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
legend("topright", paste("NORMAL",length(which(dmc_table$mean.diff>.25 & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
legend("topleft", paste("TUMOR",length(which(dmc_table$mean.diff<(-.25) & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
dev.off()
##########
saveRDS(x,"heatmap_diff_TUMOR_VS_NORMAL_x_clust.rds")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

meth.norm.sig$mean.diff>0

dim( annotation.table[rownames(annotation.table) %in% rownames(meth.norm.sig),]  )


anno_c1= annotation.table[rownames(annotation.table) %in% c1,c(1,2,2)]
anno_c1[,2] = anno_c1[,2]-150
