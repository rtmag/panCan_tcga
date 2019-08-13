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
system("mkdir TCGA_RNA_SEQ_counts/"
for(file in dir("./TCGA_RNA_SEQ_counts/",pattern = "counts.gz", recursive = T)){
    TCGAbiolinks::move(file,basename(file))
}
