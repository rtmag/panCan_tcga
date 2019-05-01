library(umap)
mval.sig = readRDS("mval_TUMORonly_tumor_vs_normal_FDR5p.rds")

mval.umap = umap(mval.sig)
