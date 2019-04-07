setwd("/root/TCGA/Rnbeads")

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
  
  # Mutation matrix
  mut.file.ix <- mut.file[mut.file$Tumor_Sample_Barcode %in% mut.id.withMutation,]
  
  
  sort(table(mut.file$Variant_Classification))
  
  tail(sort(table(paste(mut.file$Hugo_Symbol,mut.file$Tumor_Sample_Barcode))))
  
  
  # CNA table
  
  # Clinical table
  
  
  
##############################

# Sample annotation
sample.annotation <- file.path("/home/rtm/vivek/navi/EPIC/rnbeads_sample_sheet.csv")
rnb.options(import.table.separator=",")

# Report directory
system("rm -fr /home/rtm/vivek/navi/EPIC/RnBeads/RnBeads_normalization")
report.dir <- file.path("/home/rtm/vivek/navi/EPIC/RnBeads/RnBeads_normalization")

# Vanilla parameters
rnb.options(identifiers.column="Sample_ID")

# Multiprocess
num.cores <- 20
parallel.setup(num.cores)

data.source <- c(idat.dir, sample.annotation)
result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)
rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="enmix.oob")

save.rnb.set(rnb.set.norm,path="/home/rtm/vivek/navi/EPIC/RnBeads/RnBeads_normalization/rnb.set.norm.RData")
