setwd(snakemake@params[["wd"]])

library(Biostrings)

fasta<-readAAStringSet(unlist(snakemake@input))
fasta<-fasta[!duplicated(as.character(fasta))]

names(fasta)=paste0("pg|",names(fasta),"|")

writeXStringSet(fasta,unlist(snakemake@output))