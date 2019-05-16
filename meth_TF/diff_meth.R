# Filters, differential CPG between tumor and normal
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
library(graphics)

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
#
table(abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05)
library(gplots)
library(factoextra)
library(RColorBrewer)

meth.norm = meth(rnb.set.norm,row.names=T)
meth.norm.sig=meth.norm[abs(dmc_table$mean.diff)>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]
dim(meth.norm.sig)
meth.norm.sig= meth.norm.sig[complete.cases(meth.norm.sig),]
dim(meth.norm.sig)
meth.norm.sig= meth.norm.sig[apply(meth.norm.sig,1,sd)>0,]
dim(meth.norm.sig)

track= TUMOR 
track[track=="TUMOR"]=1
track[track=="NORMAL"]=2
track=as.numeric(track)
colores=c("red","white")
clab=as.character(colores[track])

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

png("heatmap_diff_TUMOR_VS_NORMAL_FDR5p.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TUMOR","NORMAL"),fill=c("red","white"), border=T, bty="n" )
dev.off()

saveRDS(x,"heatmap_diff_TUMOR_VS_NORMAL_x_clust.rds")

pdf("volcano_diff_TUMOR_VS_NORMAL.pdf")
smoothScatter(dmc_table$mean.diff,-log10(dmc_table$diffmeth.p.adj.fdr),
              ylab=expression('-Log'[10]*' Q-values'), xlab="Differential Beta-score" )
  
abline(v=-.25,lty = 2,col="grey")
abline(v=.25,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
legend("topright", paste("NORMAL",length(which(dmc_table$mean.diff>.25 & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
legend("topleft", paste("TUMOR",length(which(dmc_table$mean.diff<(-.25) & dmc_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
dev.off()


hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=6 )

track=as.numeric(groups)
colores=c("red","blue","green","orange","purple","black")
clab=(colores[track])
png("heatmap_diff_TUMOR_VS_NORMAL_FDR5p_6K.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
 heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
dev.off()

######################################################################################################################
library(ComplexHeatmap)
mut_sig = readRDS("/root/TCGA/1_COAD_test/mut_matrix_10p_filtered.rds")
rownames(mut_sig) = rownames(rnb.set.norm@pheno)[rnb.set.norm@pheno$Tumor=="TUMOR"]

x = mut_sig
x[x==1] = "MUT"
x[x==0] = ""

normal_prot = matrix("NORMAL", sum(rnb.set.norm@pheno$Tumor=="NORMAL"), dim(x)[2])
colnames(normal_prot) = colnames(x)
rownames(normal_prot) = rownames(rnb.set.norm@pheno)[rnb.set.norm@pheno$Tumor=="NORMAL"]
x = rbind(x ,normal_prot )

hp = readRDS("heatmap_diff_TUMOR_VS_NORMAL_x_clust.rds")

library(ComplexHeatmap)

col = c("MUT" = "red", "NORMAL"="white")

alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # bug red
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
            gp = gpar(fill = col["MUT"], col = NA))
    },
    # white
    NORMAL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
            gp = gpar(fill = col["NORMAL"], col = NA))
    }
)

heatmap_legend_param = list(title = "Alterations", at = c( "MUT","NORMAL"), 
        labels = c("Mutation, Indels","NORMAL"))

pdf("oncoprint.pdf")
oncoPrint(t(x),
    alter_fun = alter_fun, col = col, column_order=hp$colInd,
              row_names_gp = gpar(fontsize = 2),
     heatmap_legend_param = heatmap_legend_param)
dev.off()

# only tp53, kras, braf
three_genes = x[colnames(x) %in% c("TP53","KRAS","BRAF")]
pdf("oncoprint_tp53_braf_kras.pdf")
oncoPrint(t(three_genes),
    alter_fun = alter_fun, col = col, column_order=hp$colInd,
     heatmap_legend_param = heatmap_legend_param)
dev.off()
##################################################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(RTCGA.clinical)

clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")

coad_clinical = data.frame(times = clinical$days_to_last_follow_up,
                           bcr_patient_barcode = clinical$bcr_patient_barcode,
                           patient.vital_status = clinical$vital_status)

coad_clinical$patient.vital_status = as.character(coad_clinical$patient.vital_status)

coad_clinical$times[coad_clinical$patient.vital_status == "dead"] <- clinical$days_to_death[coad_clinical$patient.vital_status == "dead"]
# alive=0 and dead=1
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="alive"] = 0
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="dead"] = 1

groups=cutree( hc, k=2 )
meth.k.id <- data.frame( do.call( rbind, strsplit( names(groups), '-' ) ) )
meth.k.id <- paste(meth.k.id[,1],meth.k.id[,2],meth.k.id[,3],sep="-")

coad_clinical.meth <- coad_clinical[coad_clinical$bcr_patient_barcode %in% meth.k.id,]
coad_clinical.meth <- data.frame(coad_clinical.meth,Methylation.group=groups[match(coad_clinical.meth$bcr_patient_barcode,meth.k.id)])
coad_clinical.meth$Methylation.group <- as.factor(coad_clinical.meth$Methylation.group)
coad_clinical.meth$patient.vital_status <- as.numeric(coad_clinical.meth$patient.vital_status)

pdf("COAD_methylation_survival_2G.pdf")
kmTCGA(coad_clinical.meth, explanatory.names = "Methylation.group",  pval = TRUE, risk.table=FALSE,
      palette = c("red","blue"))
dev.off()

##################################################################################################################
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

hi_normal = meth.norm[dmc_table$mean.diff>.25 & dmc_table$diffmeth.p.adj.fdr<0.05,]
hi_normal= hi_normal[complete.cases(hi_normal),]
hi_normal= hi_normal[apply(hi_normal,1,sd)>0,]

hi_tumor = meth.norm[dmc_table$mean.diff<(-.25) & dmc_table$diffmeth.p.adj.fdr<0.05,]
hi_tumor= hi_tumor[complete.cases(hi_tumor),]
hi_tumor= hi_tumor[apply(hi_tumor,1,sd)>0,]

hi_tumor_anno = annotation.table[rownames(annotation.table) %in% rownames(hi_tumor),c(1,2,2)]
hi_tumor_anno[,2] = hi_tumor_anno[,2]-50
hi_tumor_anno[,3] = hi_tumor_anno[,3]+50
write.table(hi_tumor_anno, "hi_tumor_cpg.bed",sep="\t",quote=F,row.names=F,col.names=F)

hi_normal_anno = annotation.table[rownames(annotation.table) %in% rownames(hi_normal),c(1,2,2)]
hi_normal_anno[,2] = hi_normal_anno[,2]-50
hi_normal_anno[,3] = hi_normal_anno[,3]+50
write.table(hi_normal_anno, "hi_normal_cpg.bed.bed",sep="\t",quote=F,row.names=F,col.names=F)

