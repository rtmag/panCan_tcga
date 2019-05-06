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
library(betareg)
source("https://raw.githubusercontent.com/rtmag/betareg-stepwise-AIC/master/stepwise_aic_betareg.R")

# Fit the full model 
full.model <- glm(cpg ~ . ,family = quasibinomial,data=cpg)
full.model <- glm(cpg ~ . ,family = binomial(link = "logit"),data=cpg)
mod <- betareg(cpg ~ . , data = cpg, link = "logit")
plot(mod)


x = "mod <- betareg(cpg ~ . , data = cpg, link = 'logit')"
eval(parse(text=x))

gene_names = colnames(mut_sig)

x = paste("mod <- betareg(cpg ~",
      paste(colnames(mut_sig),collapse = '+'),
      ", data = cpg, link = 'logit')")
eval(parse(text=x))
AIC(mod)

x = paste("mod <- betareg(cpg ~",
      paste(colnames(mut_sig)[-171],collapse = '+'),
      ", data = cpg, link = 'logit')")
eval(parse(text=x))


current_covar = colnames(mut_sig)
model <- betareg(cpg ~ . , data = cpg, link = "logit")
sumMod = summary(model)

continue = TRUE
while(continue == TRUE){
  which.max(sumMod$coefficients$mean[-1,4])
  }


eval(parse(text="xtest = 1"))

boxplot(cpg ~ APC,data=cpg)


pval = matrix(1,nrow=1000,ncol=dim(mut_sig)[2])
rownames(pvals) = rownames(beta.sig)[1:1000]
colnames(pvals) = colnames(mut_sig)

cpg = beta.sig[1,]
cpg = data.frame(cpg = cpg, mut_sig)
for(i in 1:1000){
  cpg$cpg = beta.sig[i,]
#  model <- glm(cpg ~ . ,family = quasibinomial,data=cpg)
  model <- betareg(cpg ~ . , data = cpg, link = "logit")
  sumMod = summary(model)
#  pval[i,] = sumMod$coefficients[-1,4]
   pval[i,] = sumMod$coefficients$mean[-1,4]
}
#FDR correction
for(i in 1:dim(mut_sig)[2]){
 pval[,i] <- p.adjust(pval[,i], "BH")
}

ix = which(pval<0.05,arr.ind = TRUE)
cpg$cpg = beta.sig[which(pvals<0.05,arr.ind = TRUE)[1],]
colnames(pval)[which(pval<0.05,arr.ind = TRUE)[2]]
cxcpg = data.frame(beta=beta.sig[ix[1,1],], mut=mut_sig[,ix[1,2]])

mut_ix = 20
sub_ix = ix[ix[,2]==mut_ix,]
par(mfrow=c(3,3))
for(i in 1:9){
  cxcpg = data.frame(beta=beta.sig[sub_ix[i,1],], mut=mut_sig[,sub_ix[i,2]])
  boxplot(beta ~ mut ,cxcpg)
}

#############################################################################################
library(betareg)
beta.sig = readRDS("beta_TUMORonly_tumor_vs_normal_FDR5p.rds")
mut_sig = readRDS("~/CSI/TCGA_constrained_inference/mut_matrix_10p_filtered.rds")

pval = matrix(1,nrow=dim(beta.sig)[1],ncol=dim(mut_sig)[2])
rownames(pval) = rownames(beta.sig)
colnames(pval) = colnames(mut_sig)

cpg = beta.sig[1,]
cpg = data.frame(cpg = cpg, mut_sig)
for(i in 1:dim(beta.sig)[1]){
  cpg$cpg = beta.sig[i,]
  current_covar = colnames(mut_sig)
  x = paste("model <- betareg(cpg ~",
      paste(current_covar,collapse = '+'),
      ", data = cpg, link = 'logit')")
    
  eval(parse(text=x))
  prevAIC = AIC(model)
  
  continue = TRUE
  while(continue == TRUE){
    
      sumMod = summary(model)
      current_covar = current_covar[-which.max(sumMod$coefficients$mean[-1,4])]
      
      x = paste("model_updated <- betareg(cpg ~",
          paste(current_covar,collapse = '+'),
          ", data = cpg, link = 'logit')")
    
      eval(parse(text=x))
      updatedAIC = AIC(model_updated)
       
      if( updatedAIC <= prevAIC ){
            prevAIC = updatedAIC
            model = model_updated
        }
    
      if( updatedAIC > prevAIC ){
            continue = FALSE            
        }
  }
  
  final_pval = sumMod$coefficients$mean[-1,4]
  for( j in 1:length(final_pval)){ 
      pval[i,which(colnames(pval) == gsub("1$","",names(final_pval[j]),perl=TRUE))] = final_pval[j]
  }
  print(paste("analysing CpG:",i,"out of:",dim(beta.sig)[1]))
}

saveRDS(pval,"pval.rds")

fdr = pval
#FDR correction
for(i in 1:dim(mut_sig)[2]){
 fdr[,i] <- p.adjust(pval[,i], "BH")
}

saveRDS(fdr,"fdr.rds")

refine_test=sumMod$coefficients$mean[-1,4]
final_pval

refine_test[names(refine_test) %in% names(final_pval)]
pval_comp=data.frame(FullModel_pvals=refine_test[names(refine_test) %in% names(final_pval)],RefinedModel_pvals=final_pval)

write.csv(pval_comp,"pval_comparison_full_vs_reduced_model.csv")

#############################################################################################
# no refined
library(betareg)
beta.sig = readRDS("beta_TUMORonly_tumor_vs_normal_FDR5p.rds")
mut_sig = readRDS("~/CSI/TCGA_constrained_inference/mut_matrix_10p_filtered.rds")

pval = matrix(1,nrow=dim(beta.sig)[1],ncol=dim(mut_sig)[2])
rownames(pval) = rownames(beta.sig)
colnames(pval) = colnames(mut_sig)

cpg = beta.sig[1,]
cpg = data.frame(cpg = cpg, mut_sig)
for(i in 1:dim(beta.sig)[1]){
  cpg$cpg = beta.sig[i,]
  model <- betareg(cpg ~ ., data = cpg, link = 'logit')
  sumMod = summary(model)
  final_pval = sumMod$coefficients$mean[-1,4]
  pval[i,] = final_pval
  print(paste("analysing CpG:",i,"out of:",dim(beta.sig)[1]))
}

saveRDS(pval,"pval.rds")

fdr = pval
#FDR correction
for(i in 1:dim(mut_sig)[2]){
 fdr[,i] <- p.adjust(pval[,i], "BH")
}
saveRDS(fdr,"fdr.rds")
#############################################################################################
colnames(fdr)[as.numeric(names(sort(table((ix[,2])))))]
#############################################################################################
mx=colnames(mut_sig) %in% c("KRAS","BRAF","TP53","PCDH10","LAMA1","PLXNA4","FAT2","FBXW7","PCDH17","FRAS1")
x = mut_sig[,mx]
x2 = matrix(as.numeric(unlist(x)),ncol=10)
rownames(x2) = rownames(mut_sig)
colnames(x2) = colnames(x)
heatmap.2(t(x2),trace="none")
