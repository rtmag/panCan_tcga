library(gplots)
library(factoextra)
library(RColorBrewer)

mut_sig = readRDS("~/CSI/TCGA_constrained_inference/mut_matrix_10p_filtered.rds")

x = matrix(as.numeric(unlist(mut_sig)),ncol=dim(mut_sig)[2])
colnames(x) = colnames(mut_sig)
rownames(x) = rownames(mut_sig)

heatmap(x)
library(pheatmap)
colors <- (colorRampPalette( (brewer.pal(9, "RdYlBu")) )(200))


distance_binary=dist(x, method="binary")
sampleDistMatrix <- as.matrix(distance_binary)

 pheatmap(sampleDistMatrix,
         clustering_distance_rows=distance_binary,
         clustering_distance_cols=distance_binary,
         col=colors)
         
library(ade4)
d <- dist.binary(mut_sig, method = 1, diag = FALSE, upper = FALSE) #method 1 is Jaccard index (1901) S3 coefficient of Gower & Legendre
hc <- hclust(d)               # apply hierarchical clustering 
plot(hc)
