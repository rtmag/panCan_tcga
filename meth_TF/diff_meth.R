# Filters, differential CPG between tumor and normal
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
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

meth.norm = meth(rnb.set.norm)
meth.norm.sig=meth.norm[abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]

track= TUMOR 
track[track=="TUMOR"]=1
track[track=="NORMAL"]=2
track=as.numeric(track)
colores=c("red","white")
clab=as.character(colores[track])

png("heatmap_diff_TUMOR_VS_NORMAL_FDR5p.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TUMOR","NORMAL"),fill=c("red","white"), border=T, bty="n" )
dev.off()
##########
