more COAD_hi_normal_cpg.anno|cut -f16|sort|uniq > COAD_hi_normal_cpg.genes
more COAD_hi_tumor_cpg.anno|cut -f16|sort|uniq > COAD_hi_tumor_cpg.genes
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
library(DESeq2)

exp <- readRDS("~/CSI/TCGA/TF_meth_TCGA/COAD_RNA_SEQ.rds")
exp <- exp[rowSums(exp)>0,]

design<-data.frame(group=1:dim(exp)[2])
dds <- DESeqDataSetFromMatrix(
       countData = exp,
       colData = design,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)


# Patient Link RNA - METH
patient.anno <- read.csv("COAD_TP53_mutation_info_withNormal.csv",stringsAsFactors=F)

patient_parced <- gsub("\\-...\\-....\\-..$","",patient.anno[,1],perl=TRUE)
rna_parced <- gsub("\\-...\\-....\\-..$","",colnames(exp),perl=TRUE)
patient.anno = patient.anno[patient_parced %in% rna_parced,]
patient_parced <- gsub("\\-...\\-....\\-..$","",patient.anno[,1],perl=TRUE)

ix = match(patient_parced, rna_parced)

exp_meth = exp[,ix]
exp_norm <- exp_meth[,patient.anno$Variant_Classification=="NORMAL"]
exp_tumor <- exp_meth[,patient.anno$Variant_Classification!="NORMAL"]

################################################################################################
normal <- read.table("COAD_hi_normal_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
tumor <- read.table("COAD_hi_tumor_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
normal <- normal[,1]
tumor <- tumor[,1]


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
