mval.sig = readRDS("mval_TUMORonly_tumor_vs_normal_FDR5p.rds")
beta.sig = readRDS("beta_TUMORonly_tumor_vs_normal_FDR5p.rds")
mut_sig = readRDS("~/CSI/TCGA_constrained_inference/mut_matrix_10p_filtered.rds")


for(i in 1:dim(mut_sig)[2] ){
  mut_sig[,i] = as.factor(mut_sig[,i])
}

rownames(mut_sig) = colnames(mval.sig)


cpg = beta.sig[1,]
cpg = data.frame(cpg = cpg, mut_sig)

library(MASS)
# Fit the full model 
full.model <- glm(cpg ~ . ,family = binomial,data=cpg)


boxplot(cpg ~ APC,data=cpg)

pvals = matrix(1,nrow=1000,ncol=dim(mut_sig)[2])
rownames(pvals) = rownames(beta.sig)[1:1000]
colnames(pvals) = colnames(mut_sig)

cpg = beta.sig[1,]
cpg = data.frame(cpg = cpg, mut_sig)
for(i in 1:1000){
  cpg$cpg = beta.sig[i,]
  model <- glm(cpg ~ . ,family = binomial,data=cpg)
  sumMod = summary(model)
  pvals[i,] = sumMod$coefficients[-1,4]
}


