setwd("/root/TCGA/Rnbeads")
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)

## preprocessing
#idat files
idat.dir <- file.path("/root/TCGA/tcgaBiolink/")

#### Link list processing ####
# master link list
link.list <- read.table("/root/TCGA/tcgaBiolink/idat_filename_case.txt",header=T,sep="\t")
# remove double associated methylation profiles
link.list <- link.list[grep(",",link.list$cases,invert=TRUE),]
# get project names
projects <- sort(unique(gsub("TCGA\\-","",link.list$project,perl=TRUE)))
# Make sure it matches the panCan data
panCan.dir <- list.dirs(path = "/root/TCGA/panCancer_2018", full.names = TRUE, recursive = FALSE)
for(i in 1:length(projects)){  if(grep(tolower(projects[i]),panCan.dir) > 0){print(paste(projects[i],": OK"))}  }

# Get all genes that can have mutations across all projects
master.gene.list <- read.table(pipe("cat /root/TCGA/panCancer_2018/*/data_mutations_mskcc.txt|tail -n +2|cut -f1|sort|uniq"),sep="\n")
master.gene.list <- as.character(master.gene.list[,1])
##############################

#### Project-based processing ####
for(i in 1:length(projects)){
  #project specific link.list
  print("Parsing",projects[i],"Files")
  
  i.link.list = link.list[grep(projects[i],link.list$project),]
  
  #### subset samples with mutation information ####
  mut.file <- paste(panCan.dir[i],"/data_mutations_mskcc.txt",sep="")
  mut.file <- paste("cut -f1,10,17,40 ",mut.file)
  mut.file <- read.table(pipe(mut.file),sep="\t",header=T,quote="")
  
  meth.id.original <- unique(as.character(i.link.list$cases))
  meth.id <- data.frame( do.call( rbind, strsplit( meth.id.original, '-' ) ) )
  sample.type <- gsub("\\w$","",meth.id[,4],perl=TRUE)
  patient.id <- paste(meth.id[,1],meth.id[,2],meth.id[,3],sep="-")
  meth.id <- paste(meth.id[,1],meth.id[,2],meth.id[,3],sample.type,sep="-")
  
  meth.id.withMutation <- meth.id.original[meth.id %in% mut.file$Tumor_Sample_Barcode]
  mut.id.withMutation <- meth.id[meth.id %in% mut.file$Tumor_Sample_Barcode]
  patient.id.withMutation <- patient.id[meth.id %in% mut.file$Tumor_Sample_Barcode]
  meth.id.normal <- meth.id.original[as.numeric(sample.type)>9 & as.numeric(sample.type)<20]
  
  # create project directory
  system(paste("mkdir",projects[i]))
  
  # Original mutation matrix
  mut.file.ix <- mut.file[mut.file$Tumor_Sample_Barcode %in% mut.id.withMutation,]
  write.csv(mut.file.ix,paste(projects[i],"/",projects[i],"_mutation_original_table.csv",sep=""),row.names=F)
  
  # Mutation Matrix
  mut.file.ix <- mut.file.ix[mut.file.ix$Variant_Classification %in% c("Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation",
                                            "Nonstop_Mutation","In_Frame_Del","Frame_Shift_Ins","Nonstop_Mutation"),]
  
  mut.mat <- matrix(0, nrow=length(meth.id.withMutation), ncol=length(master.gene.list) )
  rownames(mut.mat) <- meth.id.withMutation
  colnames(mut.mat) <- master.gene.list
  mut.mat.ix <- unique(mut.file.ix[,c(1,3)])  
  for(ix in 1:dim(mut.mat.ix)[1]){ mut.mat[ mut.id.withMutation %in% mut.mat.ix[ix,2] , master.gene.list %in% mut.mat.ix[ix,1] ] = 1 }
  write.csv(mut.mat,paste(projects[i],"/",projects[i],"_mutation_matrix.csv",sep=""),row.names=F)
  if(length(meth.id.normal)>0){ tmp <- matrix(5, nrow=length(meth.id.normal), ncol=length(master.gene.list) )
                                rownames(tmp) <- colnames(meth.id.normal); 
                                colnames(tmp) <- colnames(master.gene.list); 
                               mut.mat <- rbind( mut.mat, tmp) }
  write.csv(mut.mat,paste(projects[i],"/",projects[i],"_mutation_matrix_withNormal.csv",sep=""),row.names=F)
  
  # TP53 information matrix
  mut.file.p53 <- mut.file.ix[mut.file.ix$Hugo_Symbol=="TP53",]
  mut.file.p53 <- unique(mut.file.p53)
  
  mut.file.p53 <- cbind(aggregate(as.character(Variant_Classification) ~ 
                     as.character(Tumor_Sample_Barcode), data = mut.file.p53,paste, collapse = "||"),
           aggregate(as.character(HGVSp_Short) ~ 
            as.character(Tumor_Sample_Barcode), data = mut.file.p53,paste, collapse = "||") )
  
  mut.file.p53 <- data.frame(Tumor_Sample_Barcode=mut.file.p53[,1], Variant_Classification=mut.file.p53[,2],HGVSp_Short=mut.file.p53[,4] )
  mut.file.p53 <- data.frame(meth.id.withMutation,mut.file.p53[match(mut.id.withMutation, as.character(mut.file.p53[,1]) ),])
  mut.file.p53 <- mut.file.p53[,c(1,3,4)]
  mut.file.p53[,2] <- as.character(mut.file.p53[,2])
  mut.file.p53[,2][is.na(mut.file.p53[,2])] <- "WT"
  mut.file.p53[,3] <- as.character(mut.file.p53[,3])
  mut.file.p53[,3][is.na(mut.file.p53[,3])] <- "WT"
  write.csv(mut.file.p53,paste(projects[i],"/",projects[i],"_TP53_mutation_info.csv",sep=""),row.names=F)
  if(length(meth.id.normal)>0){ tmp <- cbind(meth.id.normal,"NORMAL","NORMAL");
                                colnames(tmp) <- colnames(mut.file.p53); 
                               mut.file.p53 <- rbind( mut.file.p53, tmp) }
  rownames(mut.file.p53) <- NULL
  write.csv(mut.file.p53,paste(projects[i],"/",projects[i],"_TP53_mutation_info_withNormal.csv",sep=""),row.names=F)
  
  # CNA table
  CNA <- paste(panCan.dir[i],"/data_CNA.txt",sep="")
  CNA <- read.table(CNA,sep="\t",header=TRUE,quote="")
  CNA <- CNA[!duplicated(CNA$Hugo_Symbol), ]
  rownames(CNA) <- CNA$Hugo_Symbol
  CNA <- CNA[,3:dim(CNA)[2]]
  CNA <- t(CNA)
  rownames(CNA) <- gsub("\\.","\\-",rownames(CNA),perl=TRUE)
  CNA <- CNA[rownames(CNA) %in% mut.id.withMutation,]
  write.csv(CNA,paste(projects[i],"/",projects[i],"_CNA_matrix.csv",sep=""))

  # Clinical table
  clinical <- paste(panCan.dir[i],"/data_clinical_patient.txt",sep="")
  clinical <- paste("tail -n+5",clinical,"|cut -f 1,5,6,7,8,9,10,11,26,27,28,29,30")
  clinical <- read.table(pipe(clinical),sep="\t",header=T,quote="")
  clinical.patient.id <- data.frame( do.call( rbind, strsplit( mut.id.withMutation, '-' ) ) )
  clinical.patient.id <- paste(clinical.patient.id[,1],clinical.patient.id[,2],clinical.patient.id[,3],sep="-")
  clinical <- clinical[match(clinical.patient.id, as.character(clinical$PATIENT_ID) ),]
  clinical$PATIENT_ID <- meth.id.withMutation
  write.csv(clinical,paste(projects[i],"/",projects[i],"_clinical_info.csv",sep=""),row.names=F)

  ##############################
  print("Starting 450K processing of",projects[i])


  # Sample annotation creation
  master.meth.id <- c(meth.id.withMutation,meth.id.normal)
  sentrix <- as.character(i.link.list$file_name)
  sentrix <- gsub("_Red.idat","",sentrix)
  sentrix <- gsub("_Grn.idat","",sentrix)
  sentrix <- data.frame( do.call( rbind, strsplit( sentrix, '_' ) ) )
  sentrix <- unique(cbind(i.link.list[,1],sentrix))
  sentrix <- sentrix[match(master.meth.id,sentrix[,1]),]
  sample.annotation <- data.frame(Sample_ID=as.character(sentrix[,1]),
                                  Sentrix_ID=as.character(sentrix[,2]),
                                  Sentrix_Position=as.character(sentrix[,3]) )
  write.csv(sample.annotation,paste(projects[i],"/",projects[i],"_sample_annotation.csv",sep=""),row.names=F)

  # File path
  sample.annotation <- file.path(paste(projects[i],"/",projects[i],"_sample_annotation.csv",sep=""))
  rnb.options(import.table.separator=",")

  # Report directory
  command <- paste("rm -fr ",projects[i],"/","RnBeads_normalization/",sep="")
  system(command)
  report.dir <- file.path(paste(projects[i],"/","RnBeads_normalization/",sep=""))

  # Vanilla parameters
  rnb.options(identifiers.column="Sample_ID")

  # Multiprocess
  num.cores <- 20
  parallel.setup(num.cores)
  
  #idat files
  idat.dir <- file.path("/root/TCGA/tcgaBiolink")

  data.source <- c(idat.dir, sample.annotation)
  result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)
  rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="methylumi.noob")

  save.rnb.set(rnb.set.norm,path=paste(projects[i],"/","RnBeads_normalization/rnb.set.norm_withNormal.RData",sep=""))
  
  meth.norm<-meth(rnb.set.norm)
  colnames(meth.norm) = as.character(rnb.set.norm@pheno[,1])
  rownames(meth.norm) = rownames(rnb.set.norm@sites)
  saveRDS(meth.norm, paste(projects[i],"/","RnBeads_normalization/betaVALUES_withNormal.rds",sep=""))

  mval.norm <- mval(rnb.set.norm,row.names=T)
  colnames(mval.norm) = as.character(rnb.set.norm@pheno[,1])
  rownames(mval.norm) = rownames(rnb.set.norm@sites)
  saveRDS(mval.norm, paste(projects[i],"/","RnBeads_normalization/mVALUES_withNormal.rds",sep=""))
  
  if(length(meth.id.normal)>0){
    rnb.set.norm_noNormal=remove.samples(rnb.set.norm,samples(rnb.set.filtered)[which(mut.file.p53$Variant_Classification=="NORMAL")])
    save.rnb.set(rnb.set.norm_noNormal,path=paste(projects[i],"/","RnBeads_normalization/rnb.set.norm.RData",sep=""))
    
    meth.norm<-meth(rnb.set.norm)
    colnames(meth.norm) = as.character(rnb.set.norm@pheno[,1])
    rownames(meth.norm) = rownames(rnb.set.norm@sites)
    saveRDS(meth.norm, paste(projects[i],"/","RnBeads_normalization/betaVALUES.rds",sep=""))

    mval.norm <- mval(rnb.set.norm,row.names=T)
    colnames(mval.norm) = as.character(rnb.set.norm@pheno[,1])
    rownames(mval.norm) = rownames(rnb.set.norm@sites)
    saveRDS(mval.norm, paste(projects[i],"/","RnBeads_normalization/mVALUES.rds",sep=""))
   }
    
  if(length(meth.id.normal)<1){
    cp.command <- paste("cp ",projects[i],"/","RnBeads_normalization/rnb.set.norm_withNormal.RData ",projects[i],"/",
          "RnBeads_normalization/rnb.set.norm.RData",sep="")
    system(cp.command)
    
    cp.command <- paste("cp ",projects[i],"/","RnBeads_normalization/betaVALUES_withNormal.rds ",projects[i],"/",
          "RnBeads_normalization/betaVALUES.rds",sep="")
    system(cp.command)
    
    cp.command <- paste("cp ",projects[i],"/","RnBeads_normalization/mVALUES_withNormal.rds ",projects[i],"/",
          "RnBeads_normalization/mVALUES.rds",sep="")
    system(cp.command)
    }
}


########################
