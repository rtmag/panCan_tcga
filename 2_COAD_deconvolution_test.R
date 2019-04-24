suppressMessages(library(RnBeads))
library("ChAMP")
options(bitmapType="cairo")
options(scipen=999)
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


