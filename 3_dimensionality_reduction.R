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
library(umap)

