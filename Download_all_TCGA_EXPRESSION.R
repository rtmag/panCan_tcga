
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL
for(proj in projects){
    print(proj)
    query <- GDCquery(project = proj,
                      data.category = "Raw microarray data",
                      data.type = "Raw intensities", 
                      experimental.strategy = "Methylation array", 
                      legacy = TRUE,
                      file.type = ".idat",
                      platform = "Illumina Human Methylation 450")
    match.file.cases <- getResults(query,cols=c("cases","file_name"))
    match.file.cases$project <- proj
    match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
             error = function(e) GDCdownload(query, method = "client"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
    TCGAbiolinks::move(file,basename(file))
}

query.exp.hg38 <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ",
                  barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp.rda")

