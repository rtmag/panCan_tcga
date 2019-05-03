library(umap)
mval.sig = readRDS("mval_TUMORonly_tumor_vs_normal_FDR5p.rds")

mval.sd = apply(mval.sig,1,sd)

mval.umap = umap(mval.sig)
saveRDS(mval.umap,"umap.rds")
mval.umap2 = umap(t(mval.sig[mval.sd>1.5,]))

pdf("umap.pdf")
plot(mval.umap2$layout[,1], mval.umap2$layout[,2])
dev.off()

library(Rtsne)
tsne = Rtsne(t(mval.sig))
saveRDS(tsne,"tsne.rds")
pdf("test_tsne.pdf")
plot(tsne$Y)
dev.off()

library(NbClust)
res<-NbClust(mval.sig, distance = "euclidean",method = "kmeans")
saveRDS(res,"res.rds")

KMeansSparseCluster(x, K=)

####################################################################
beta = readRDS("beta_TUMORonly_tumor_vs_normal_FDR5p.rds")

mval = readRDS("mval_TUMORonly_tumor_vs_normal_FDR5p.rds")
mval = mval[complete.cases(mval), ]

pdf("pca_parameters.pdf")
par(mfrow=c(2,2))

ir.pca <- prcomp(t(mval),
                 center = FALSE,
                 scale. = FALSE) 
sx=summary(ir.pca)
plot(ir.pca$x[,1],ir.pca$x[,2],xlab=paste("PC1:",round(sx$importance[2,1]*100,digits=1),"%"),
     ylab=paste("PC2:",round(sx$importance[2,2]*100,digits=1),"%"),pch=19)
#########
ir.pca <- prcomp(t(mval),
                 center = FALSE,
                 scale. = TRUE) 
sx=summary(ir.pca)
plot(ir.pca$x[,1],ir.pca$x[,2],xlab=paste("PC1:",round(sx$importance[2,1]*100,digits=1),"%"),
     ylab=paste("PC2:",round(sx$importance[2,2]*100,digits=1),"%"),pch=19)
#########
ir.pca <- prcomp(t(mval),
                 center = TRUE,
                 scale. = FALSE) 
sx=summary(ir.pca)
plot(ir.pca$x[,1],ir.pca$x[,2],xlab=paste("PC1:",round(sx$importance[2,1]*100,digits=1),"%"),
     ylab=paste("PC2:",round(sx$importance[2,2]*100,digits=1),"%"),pch=19)
#########
ir.pca <- prcomp(t(mval),
                 center = TRUE,
                 scale. = TRUE) 
sx=summary(ir.pca)
plot(ir.pca$x[,1],ir.pca$x[,2],xlab=paste("PC1:",round(sx$importance[2,1]*100,digits=1),"%"),
     ylab=paste("PC2:",round(sx$importance[2,2]*100,digits=1),"%"),pch=19)
#########
dev.off()

############################################################################################################
options(bitmapType="cairo")
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

source("https://raw.githubusercontent.com/rtmag/refactor/master/R/refactor.R")
mvaltxt <- data.frame(ID=rownames(mval),mval)
mvaltxt = rbind(colnames(mvaltxt),mvaltxt)

res_2 <- refactor(mvaltxt,2,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_3 <- refactor(mvaltxt,3,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_4 <- refactor(mvaltxt,4,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_5 <- refactor(mvaltxt,5,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_6 <- refactor(mvaltxt,6,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_7 <- refactor(mvaltxt,7,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
res_8 <- refactor(mvaltxt,8,t=1000,stdth=0.04,out="res_without_normals_stdth.05")
#
results= res_2
pdf("res2_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res2_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
x=heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
results= res_3
pdf("res3_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res3_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
results= res_4
pdf("res4_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res4_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
#
results= res_5
pdf("res5_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res5_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
results= res_6
pdf("res6_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res6_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
results= res_7
pdf("res7_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res7_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#
results= res_8
pdf("res8_scatterplot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2])
dev.off()

pdf("res8_heatmap.pdf")
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()
#



results_4_05 <- refactor(mvaltxt,4,t=1000,stdth=0.06,out="res_without_normals_stdth.06")

results= results_4_05
plot(results$refactor_components[,1],results$refactor_components[,2])
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)

##############################################################################################################
saveRDS(res_2,"refactor_2K_t1000.rds")

x=heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)

hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=3 )

track=as.numeric(groups)
colores=c("red","blue","green")
clab=(colores[track])

pdf("refactor_2_peasonClust_K3.pdf")
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1,ColSideColors=clab)
dev.off()

pdf("refactor_2_peasonClust_K3_scatterPLot.pdf")
plot(results$refactor_components[,1],results$refactor_components[,2],col=clab,pch=19)
dev.off()

##############################################################################################################
# SURVIVAL ANALYSIS
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(RTCGA.clinical)

clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
datatable(clinical, filter = 'top', 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
          rownames = FALSE)
##################################################################################################################################
coad_clinical = data.frame(times = clinical$days_to_last_follow_up,
                           bcr_patient_barcode = clinical$bcr_patient_barcode,
                           patient.vital_status = clinical$vital_status)

coad_clinical$patient.vital_status = as.character(coad_clinical$patient.vital_status)

coad_clinical$times[coad_clinical$patient.vital_status == "dead"] <- clinical$days_to_death[coad_clinical$patient.vital_status == "dead"]
# alive=0 and dead=1
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="alive"] = 0
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="dead"] = 1

meth.k.id <- data.frame( do.call( rbind, strsplit( names(groups), '-' ) ) )
meth.k.id <- paste(meth.k.id[,1],meth.k.id[,2],meth.k.id[,3],sep="-")

coad_clinical.meth <- coad_clinical[coad_clinical$bcr_patient_barcode %in% meth.k.id,]
coad_clinical.meth <- data.frame(coad_clinical.meth,Methylation.group=groups[match(coad_clinical.meth$bcr_patient_barcode,meth.k.id)])
coad_clinical.meth$Methylation.group <- as.factor(coad_clinical.meth$Methylation.group)
coad_clinical.meth$patient.vital_status <- as.numeric(coad_clinical.meth$patient.vital_status)

pdf("COAD_methylation_survival.pdf")
kmTCGA(coad_clinical.meth, explanatory.names = "Methylation.group",  pval = TRUE, risk.table=FALSE)
dev.off()

