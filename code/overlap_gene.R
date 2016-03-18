library(IRanges)
library("data.table")
chroms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
load("/home/qzzh/Renbing/hg38gene/hg38gene.Rdata")
for (i in chroms){
tmp1 <- hg38gene[hg38gene$chrom == i,]
promoterRegion <- IRanges(start=tmp1$promoterstart1000,end=tmp1$promoterend1000)
I <- paste("/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/",i,".Rdata",sep="")
load(I)
peakRegion <- IRanges(start=output$start,end=output$end)
overlap <- findOverlaps(peakRegion,promoterRegion)
output$promoterRegion <- rep(0,nrow(output))
output$promoterRegion[overlap@queryHits] <- 1
O <- paste("/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/",i,".Rdata",sep="")
save(output,file=O)
message(paste(i," ","is finished ",sep=""),date())
}


