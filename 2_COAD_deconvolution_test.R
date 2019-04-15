suppressMessages(library(RnBeads))
library("ChAMP")
options(bitmapType="cairo")
options(scipen=999)
#
rnb.set.norm=load.rnb.set("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")
rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.norm)$dataset
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
save.rnb.set(rnb.set.filtered,path="rnb.set.norm.filtered.RData")
#
beta <- meth(rnb.set.filtered,row.names=TRUE)
beta <- beta[complete.cases(beta), ]

rnb.set.filtered@pheno$Sentrix_ID <- as.factor(rnb.set.filtered@pheno$Sentrix_ID)
beta.combat <- champ.runCombat(beta=beta,pd=rnb.set.filtered@pheno,batchname=c("Sentrix_ID"))
