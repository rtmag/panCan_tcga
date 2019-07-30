suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

num.cores <- 20
parallel.setup(num.cores)

normal_sample_number <- read.table(
         pipe('grep -c "NORMAL" /root/TCGA/Rnbeads/*/*_info_withNormal.csv|perl -pe "s/\\//\t/g"|perl -pe "s/\\:/\t/g"|cut -f5,7'),
         stringsAsFactors=FALSE)

normal_sample_number <- normal_sample_number[normal_sample_number[,2]>9,1]

for(i in 1:length(normal_sample_number)){
     tcga.dir <- list.dirs(path = "/root/TCGA/Rnbeads/", full.names = TRUE, recursive = FALSE)
     command <- paste0("/root/TCGA/Rnbeads/",normal_sample_number[i],"/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")
     rnb.set.norm=load.rnb.set(command)
         
     command <- paste0("/root/TCGA/Rnbeads/",normal_sample_number[i],"/",normal_sample_number[i],"_TP53_mutation_info_withNormal.csv")
     TUMOR = read.csv(command,header=TRUE)
     TUMOR = as.character(TUMOR$Variant_Classification)
     TUMOR[TUMOR!="NORMAL"] = "TUMOR"
     rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Tumor = TUMOR)
         
     rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.norm)$dataset
     rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
     rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
         
     #save filtered object
         #make dir
     command = paste0("mkdir /root/TCGA/TF_METH/",normal_sample_number[i])
     system(command)
     command = paste0("/root/TCGA/TF_METH/",normal_sample_number[i],"/rnb.set.norm.filtered.RData")
     save.rnb.set(rnb.set.filtered,path=command)

     dmc <- rnb.execute.computeDiffMeth(rnb.set.filtered, pheno.cols=c("Tumor") )
     comparison <- get.comparisons(dmc)[1]
     print(comparison)
     dmc_table <-get.table(dmc, comparison, "sites", return.data.frame=TRUE)
         
     command = paste0("/root/TCGA/TF_METH/",normal_sample_number[i],"/",normal_sample_number[i],"_dmc_table.csv")
     write.csv(dmc_table,command)
}

