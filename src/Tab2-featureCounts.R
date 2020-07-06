library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
FC <- args[1]
BAM <- args[2]
cc <- featureCounts(BAM, annot.inbuilt="hg38", nthreads=8, allowMultiOverlap=F, countMultiMappingReads=T, isPairedEnd=T)
cout <- cbind( cc$annotation$GeneID, cc$counts )
write.table(cout, FC, quote=F, row.names=F, col.names=F, sep="\t")
