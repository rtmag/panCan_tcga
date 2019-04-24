suppressMessages(library(RnBeads))
library("ChAMP")
options(bitmapType="cairo")
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
#############################################################################################################
# NAIVE TEST
source("/root/myPrograms/refactor/R/refactor.R")

rnb.set.filtered <- load.rnb.set(path="/root/TCGA/1_COAD_test/refactor/rnb.set.norm.filtered.RData.zip")
beta <- meth(rnb.set.filtered,row.names=TRUE)
beta <- beta[complete.cases(beta), ]
saveRDS(beta,"beta.filtered.rds")

betatxt <- data.frame(ID=rownames(beta),beta)
write.table(betatxt,"betaVALUES.txt",sep="\t",quote=FALSE,row.names=FALSE)

k_estimate_command = paste("python /root/myPrograms/refactor/python/estimate_k_server.py --datafile betaVALUES.txt --max_k 15")
system( k_estimate_command )

t_estimate_command = paste("python /root/myPrograms/refactor/python/estimate_t_server.py --datafile betaVALUES.txt --k 11 --numsites 5000")
system( t_estimate_command )

#  chnge for editted version because original takes long to read files
results <- refactor("betaVALUES.txt",11,t=2000,stdth=0.01,out="demo_results")

refactorCA = read.table("demo_results.out.components.txt")
refactorCpG = read.table("demo_results.out.rankedlist.txt")

plot(refactorCA[,1],refactorCA[,2])
plot(res$standard_pca[,1],res$standard_pca[,2])

beta = readRDS("beta.filtered.rds")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

all.meth.norm = beta[rownames(beta) %in% refactorCpG[1:1000,1], ]

pdf("heatmap_sigCPG.pdf")
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",cexCol=.1)
dev.off()

###################
#bad outliers 
#TCGA-D5-6536-01A-11D-1721-05
#TCGA-D5-6924-01A-11D-1926-05
# EDITED
source("https://raw.githubusercontent.com/rtmag/refactor/master/R/refactor.R")

beta = readRDS("beta.filtered.rds")
beta = beta[,!colnames(beta) %in% c("TCGA-D5-6536-01A-11D-1721-05","TCGA-D5-6924-01A-11D-1926-05")]
betatxt <- data.frame(ID=rownames(beta),beta)
write.table(betatxt,"betaVALUES_without_outliers.txt",sep="\t",quote=FALSE,row.names=FALSE)

k_estimate_command = 
paste("python /root/myPrograms/refactor/python/estimate_k_server.py --datafile betaVALUES_without_outliers.txt --max_k 15")
system( k_estimate_command )
system(" mv estimate_k_results.png estimate_k_beta_without_outliers.png ")

t_estimate_command = 
paste("python /root/myPrograms/refactor/python/estimate_t_server.py --datafile betaVALUES_without_outliers.txt --k 11 --numsites 5000")
system( t_estimate_command )
system(" mv estimate_t_results.png estimate_t_beta_without_outliers.png ")

#########
betatxt1 = rbind(colnames(betatxt),betatxt)
results <- refactor(betatxt1,11,t=2000,stdth=0.01,out="refactor_without_outliers")


############################################################################################################
######### TEST REFACTOR
# ORIGINAL
source("/root/myPrograms/refactor/R/refactor.R")
# EDITED
source("https://raw.githubusercontent.com/rtmag/refactor/master/R/refactor.R")
#########
ORIGINAL_res <- refactor("demo_datafile.txt",2,t=500,stdth=0.01,out="demo_results_fast_ori")

beta = read.table("demo_datafile.txt",sep="\t",stringsAsFactors=FALSE)
EDITED_res <- refactor(beta,2,t=500,stdth=0.01,out="demo_results_fast")


############################################################################################################
rnb.set.norm=load.rnb.set("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")

TUMOR = read.csv("/root/TCGA/Rnbeads/COAD/COAD_TP53_mutation_info_withNormal.csv",header=TRUE)
TUMOR = as.character(TUMOR$Variant_Classification)
TUMOR[TUMOR!="NORMAL"] = "TUMOR"
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Tumor = TUMOR)

rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.norm)$dataset
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
save.rnb.set(rnb.set.filtered,path="/root/TCGA/1_COAD_test/refactor/rnb.set.norm.filtered.RData")
#
beta <- meth(rnb.set.filtered,row.names=TRUE)
beta <- beta[complete.cases(beta), ]

rnb.set.norm@pheno$Sentrix_ID <- as.factor(rnb.set.norm@pheno$Sentrix_ID)
beta.combat <- champ.runCombat(beta=beta,pd=rnb.set.norm@pheno,batchname=c("Sentrix_ID"),variablename="Tumor")

myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"),variablename="Sample_Group")

meth.id <- data.frame( do.call( rbind, strsplit( as.character(rnb.set.norm@pheno$Sample_ID), '-' ) ) )
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, plateID = meth.id[,6])
beta.combat <- champ.runCombat(beta=beta,pd=rnb.set.norm@pheno,batchname="Sentrix_ID",variablename="plateID")

# test for file removal 
x = read.table(pipe("grep -v ',' idat_filename_case.txt"),sep="\t",header=T)
meth.id <- data.frame( do.call( rbind, strsplit( as.character(x$cases), '-' ) ) )

sentrix = data.frame( do.call( rbind, strsplit( as.character(x$file_name), '_' ) ) )
sentrix = as.character(sentrix[,1])


