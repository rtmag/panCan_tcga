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

