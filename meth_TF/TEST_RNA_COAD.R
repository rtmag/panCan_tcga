more COAD_hi_normal_cpg.anno|cut -f16|sort|uniq > COAD_hi_normal_cpg.genes
more COAD_hi_tumor_cpg.anno|cut -f16|sort|uniq > COAD_hi_tumor_cpg.genes
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

patient.anno <- read.csv("COAD_TP53_mutation_info_withNormal.csv",stringsAsFactors=F)
exp <- read.table("data_RNA_Seq_v2_mRNA_median_Zscores.txt",sep="\t",header=T,stringsAsFactors=F)
normal <- read.table("COAD_hi_normal_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
tumor <- read.table("COAD_hi_tumor_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
normal <- normal[,1]
tumor <- tumor[,1]

rownames(exp) = make.names(exp[,1],unique=T)
exp = exp[,3:dim(exp)[2]]

patient_parced <- gsub(".\\-...\\-....\\-..$","",patient.anno[,1],perl=TRUE)
patient_parced <- gsub("-",".",patient_parced,perl=TRUE)

ix <- match(colnames(exp), patient_parced)

exp = exp[,!is.na(ix)]
#patient.anno[ix[!is.na(ix)],]

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))
sig_vsd = exp[rownames(exp) %in% tumor,]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
 xlab="", ylab="Genes",key.title="Gene expression",cexCol=.65,cexRow=.2)
 
tumor_u=tumor[!tumor %in% normal]
normal_u=normal[!normal %in% tumor]
 
 sig_vsd_n = exp[rownames(exp) %in% normal_u,]
  sig_vsd_t = exp[rownames(exp) %in% tumor_u,]
 boxplot(rowMeans(sig_vsd_n),rowMeans(sig_vsd_t))
