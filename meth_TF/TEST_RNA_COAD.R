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

tmpgroup <- patient.anno$Variant_Classification!="NORMAL"
tmpgroup[tmpgroup=="FALSE"] = "NORMAL"
tmpgroup[tmpgroup=="TRUE"] = "TUMOR"

design<-data.frame(group=tmpgroup)
dds <- DESeqDataSetFromMatrix(
       countData = exp_meth,
       colData = design,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dds)
vsd_meth <- assay(dLRT_vsd)
saveRDS(vsd_meth,"COAD_vsd_meth.rds")
vsd_norm <- vsd_meth[,patient.anno$Variant_Classification=="NORMAL"]
vsd_tumor <- vsd_meth[,patient.anno$Variant_Classification!="NORMAL"]
################################################################################################
normal <- read.table("COAD_hi_normal_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
tumor <- read.table("COAD_hi_tumor_cpg.genes",sep="\t",header=F,stringsAsFactors=F)
normal <- normal[,1]
tumor <- tumor[,1]


colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

tumor_u=tumor[!tumor %in% normal]
normal_u=normal[!normal %in% tumor]

sig_vsd_norm_tumor = vsd_norm[rownames(vsd_norm) %in% tumor_u,]
sig_vsd_norm_tumor = sig_vsd_norm_tumor[complete.cases(sig_vsd_norm_tumor),]

sig_vsd_tumor_tumor = vsd_tumor[rownames(vsd_tumor) %in% tumor_u,]
sig_vsd_tumor_tumor = sig_vsd_tumor_tumor[complete.cases(sig_vsd_tumor_tumor),]

sig_vsd_norm_normal = vsd_norm[rownames(vsd_norm) %in% normal_u,]
sig_vsd_norm_normal = sig_vsd_norm_normal[complete.cases(sig_vsd_norm_normal),]

sig_vsd_tumor_normal = vsd_tumor[rownames(vsd_tumor) %in% normal_u,]
sig_vsd_tumor_normal = sig_vsd_tumor_normal[complete.cases(sig_vsd_tumor_normal),]

pdf("boxplot.pdf")
 boxplot(rowMeans(sig_vsd_norm_tumor),rowMeans(sig_vsd_tumor_tumor),rowMeans(sig_vsd_norm_normal),rowMeans(sig_vsd_tumor_normal),
        names=c("NORMAL","TUMOR","NORMAL","TUMOR"),col=c("#ffb3ba","#ffb3ba","#bae1ff","#bae1ff"),ylab="Variance Stabilized Log2 Normalized Read Counts")
dev.off()

pdf("boxplot_legend.pdf")
plot.new()
legend("topright",legend=c("Genes assoc Hi-meth Tumor","Genes assoc Hi-meth Normal"),fill=c("#ffb3ba","#bae1ff"), border=T, bty="n" )
dev.off()



  heatmap.2(as.matrix(sig_vsd),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
 xlab="", ylab="Genes",key.title="Gene expression",cexCol=.65,cexRow=.2)
 
