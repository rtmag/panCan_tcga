library(TCGAbiolinks)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL

for(proj in projects){
    print(proj)
    query <- GDCquery(project = proj,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  )
    match.file.cases <- getResults(query,cols=c("cases","file_name"))
    match.file.cases$project <- proj
    match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
             error = function(e) GDCdownload(query, method = "client"))
}
             
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "RNA_SEQ_filename_case.txt")
# code to move all files to local folder
system("mkdir TCGA_RNA_SEQ_counts/")
for(file in dir("./GDCdata/",pattern = "counts.gz", recursive = T)){
#    TCGAbiolinks::move(file,paste0("TCGA_RNA_SEQ_counts/",basename(file)))
    command <- paste0(file," TCGA_RNA_SEQ_counts/")
}

       
       projects[11]
       
        query <- GDCquery(project = projects[11],
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  )
    match.file.cases <- getResults(query,cols=c("cases","file_name"))
    match.file.cases$project <- proj
    match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
             error = function(e) GDCdownload(query, method = "client"))

             #putting all together
             paste -d '=' *.counts | sed 's/=\S*//g' > ../COAD_MATRIX_pre.txt
              ls -1 > ../COAD_file_names.txt
             
             x = read.table("COAD_MATRIX_pre.txt",sep="\t",row.names=1,header=FALSE)
             xname = read.table("COAD_file_names.txt",sep="\t")
             colnames(x) = xname[,1]
             x = x[1:(dim(x)[1]-5),]
             
             

source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
             
gene.df <- bitr(rownames(x), fromType = "ENSEMBL",
                        toType = c( "ENTREZID", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)
             
getBM(attributes='hgnc_symbol', 
      filters = 'ensembl_gene_id', 
      values = rownames(x), 
      mart = ensembl)
             


library(ensembldb)
library("EnsDb.Hsapiens.v86")   
edb <- EnsDb.Hsapiens.v86
tx <- genes(edb, columns=c("gene_id", "gene_name"))
tx=data.frame(gene_id=tx$gene_id, gene_name=tx$gene_name)

             
xx = gsub("\\..+","",rownames(x),perl=T)

rownames(x)=make.names(as.character(tx[ match(xx,tx[,1]),2]),unique=T)
             
 link = read.table("RNA_SEQ_filename_case.txt",sep="\t",header=T,stringsAsFactors=F)
link$file_name = gsub(".gz","",link$file_name)
colnames(x) = link[match(colnames(x), link$file_name),1]
saveRDS(x,"COAD_RNA_SEQ.rds")
