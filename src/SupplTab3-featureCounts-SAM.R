library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
BAM <- args[1]
cc <- featureCounts(BAM, annot.inbuilt="hg38", nthreads=8, allowMultiOverlap=T, countMultiMappingReads=T, reportReads="SAM")

