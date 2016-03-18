dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
file_path1="/home/qzzh/J.R.Ecker/SB_03/SB_03_unique.mtbr.cg/"
file_path2="/home/qzzh/cgDensity/SB_03/"
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
library(methyutils)
load(O)
cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
chr <- unique(cg.mtbr$chrom)
library(BSgenome.Hsapiens.UCSC.hg38)
dna.seq <- Hsapiens[[chr]]
cposition <- GetCcontextPosition(dna.seq)
win <- list(L = 500, R = 500)
cgs <- swsCalc(cposition, win)
name<-paste("variableStep chrom=",z," span=1",sep="")
CGS<-data.frame(1:length(cgs),cgs)
cgs_p<-CGS[which(cgs!=0),]
cgs_n<-cgs_p
cgs_n[,2] <- -cgs_n[,2]
Q<-paste(file_path2,"cgs_n.txt",sep="")
W<-paste(file_path2,"cgs_p.txt",sep="")
write.table(name,file=Q,append=T,quote=F,row.names=F,col.names=F)
write.table(name,file=W,append=T,quote=F,row.names=F,col.names=F)
write.table(cgs_n,file=Q,append=T,quote=F,row.names=F,col.names=F)
write.table(cgs_p,file=W,append=T,quote=F,row.names=F,col.names=F)
}






