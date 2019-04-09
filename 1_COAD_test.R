# TUMOR AND NORMAL PCA
mval = readRDS("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/mVALUES_withNormal.rds")

mval.complete = mval[complete.cases(mval), ]

ir.pca <- prcomp(t(mval.complete),
                 center = TRUE,
                 scale. = FALSE) 

sx=summary(ir.pca)

TUMOR = read.csv("/root/TCGA/Rnbeads/COAD/COAD_TP53_mutation_info_withNormal.csv",header=TRUE)
TUMOR = as.character(TUMOR$Variant_Classification)
TUMOR[TUMOR!="NORMAL"] = "TUMOR"

track= TUMOR 
track[track=="TUMOR"]=1
track[track=="NORMAL"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

pdf("450K_methylation_allProbes_mvalues_pca.pdf")
plot(ir.pca$x[,1],ir.pca$x[,2],col="",xlab=paste("PCA1:",round(sx$importance[2,1]*100,digits=1),"%"),
     ylab=paste("PCA2:",round(sx$importance[2,2]*100,digits=1),"%"))
legend("topright",legend=c("TUMOR","NORMAL"),fill=c("#ffb3ba","#baffc9"), border=T, bty="n" )
dev.off()


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

# QC
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

track=as.character(rnb.set.norm@pheno$Tumor)
track[track=="TUMOR"]=1
track[track=="NORMAL"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

meth.norm = readRDS("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/mVALUES_withNormal.rds")
meth.norm.sig=meth.norm[which(dmc_table$diffmeth.p.adj.fdr<0.05),]
meth.norm.sig = meth.norm.sig[complete.cases(meth.norm.sig),]
saveRDS(meth.norm.sig,"beta_tumor_vs_normal_FDR5p.rds")

mval.sig=mval[which(dmc_table$diffmeth.p.adj.fdr<0.05),] # 228,166 cpg (total 485,577 cpg) - 46.98% cpg retained
mval.sig = mval.sig[complete.cases(mval.sig),]
saveRDS(mval.sig,"mval_tumor_vs_normal_FDR5p.rds")

png("heatmap_TUMOR_VS_NORMAL_FDR5p.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TUMOR","NORMAL"),fill=c("#ffb3ba","#baffc9"), border=T, bty="n" )
dev.off()
##################################

# Filter mutations that occur in at least 10% of the population
mut = read.csv("/root/TCGA/Rnbeads/COAD/COAD_mutation_matrix.csv",header=TRUE)
mut.csum = colSums(mut)
mut.10p = sort(mut.csum[ (mut.csum / dim(mut)[1]) > .10 ])
mut_sig = mut[ , colnames(mut) %in% names(mut.10p)]
