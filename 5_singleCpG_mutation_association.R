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
gene_sorted = sort(table((ix[,2])))
names_genes=colnames(fdr)[as.numeric(names(sort(table((ix[,2])))))]

names(gene_sorted) = names_genes
names(tail(gene_sorted,n=25))

associated_genes = gene_sorted[names(gene_sorted) %in% names(tail(gene_sorted,n=25))]
pdf("top25_associated_genes.pdf",width=10,height=4)
barplot(associated_genes,las=2,cex.names=1,ylab="Number of cytosines associated with a gene",cex.axis=.5)
dev.off()
#############################################################################################
mx=colnames(mut_sig) %in% names(tail(gene_sorted,n=25))
x = mut_sig[,mx]
x[x==1] = "MUT"
x[x==0] = ""

library(ComplexHeatmap)
col = c("MUT" = "red")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # bug red
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
            gp = gpar(fill = col["MUT"], col = NA))
    }
)

heatmap_legend_param = list(title = "Alterations", at = c( "MUT"), 
        labels = c("Mutation"))
pdf("oncoprint.pdf")
oncoPrint(t(x),
    alter_fun = alter_fun, col = col, 
     heatmap_legend_param = heatmap_legend_param)
dev.off()
