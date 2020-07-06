library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
R1 <- args[1]
R2 <- args[2]
if(is.na(R2) || R2=="NA") R2<-NULL
BAM <- args[3]

idx.name <- "Idx-hg38-full"

align( idx.name, R1, R2 , output_file=BAM, nthreads=8, annot.inbuilt = "hg38", useAnnotation=T)

