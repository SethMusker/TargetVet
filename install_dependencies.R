args <- commandArgs(trailingOnly = TRUE)
if(args[1]!=""){
.libPaths(args[1])
cat("Packages will be installed in",args[1],"\n")
} else{
    cat("Packages will be installed in",.libPaths()[1],"\n")
}

packs<-suppressMessages(scan("cran_deps.txt",what="char"))
cat("installing CRAN packages:\n",packs,"\n")
for(i in packs){
    # cat(i,"\n")
    if(!library(i,character.only = T,logical.return = T)) install.packages(i,repos="https://cloud.r-project.org")
}

packs<-suppressMessages(scan("BiocManager_deps.txt",what="char"))
cat("installing BiocManager packages:\n",packs,"\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
BiocManager::install("Biostrings",ask=FALSE)
BiocManager::install("Rsamtools",ask=FALSE)
