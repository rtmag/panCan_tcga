library(umap)
mval.sig = readRDS("mval_TUMORonly_tumor_vs_normal_FDR5p.rds")

mval.sd = apply(mval.sig,1,sd)

mval.umap = umap(mval.sig)
mval.umap2 = umap(t(mval.sig[mval.sd>1.5,]))

pdf("umap.pdf")
plot(mval.umap2$layout[,1], mval.umap2$layout[,2])
dev.off()

library(Rtsne)
tsne = Rtsne(t(mval.sig))
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

results= results_05



results= results_05
plot(results$refactor_components[,1],results$refactor_components[,2])
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)

results_4_05 <- refactor(mvaltxt,4,t=1000,stdth=0.06,out="res_without_normals_stdth.06")

results= results_4_05
plot(results$refactor_components[,1],results$refactor_components[,2])
all.meth.norm = beta[rownames(beta) %in% results$RankedProbeNames[1:1000], ]
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)

