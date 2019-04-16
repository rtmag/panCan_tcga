suppressMessages(library(RnBeads))
library("ChAMP")
options(bitmapType="cairo")
options(scipen=999)
#
rnb.set.norm=load.rnb.set("/root/TCGA/Rnbeads/COAD/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")

TUMOR = read.csv("/root/TCGA/Rnbeads/COAD/COAD_TP53_mutation_info_withNormal.csv",header=TRUE)
TUMOR = as.character(TUMOR$Variant_Classification)
TUMOR[TUMOR!="NORMAL"] = "TUMOR"
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Tumor = TUMOR)

rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.norm)$dataset
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
save.rnb.set(rnb.set.filtered,path="rnb.set.norm.filtered.RData")
#
beta <- meth(rnb.set.filtered,row.names=TRUE)
beta <- beta[complete.cases(beta), ]

rnb.set.norm@pheno$Sentrix_ID <- as.factor(rnb.set.norm@pheno$Sentrix_ID)
beta.combat <- champ.runCombat(beta=beta,pd=rnb.set.norm@pheno,batchname=c("Sentrix_ID"),variablename="Tumor")

myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"),variablename="Sample_Group")
