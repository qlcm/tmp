dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
file_path1="/home/qzzh/J.R.Ecker/mtbr/SB_03/SB_03_unique.mtbr.cg/"
file_path2="/home/qzzh/cgDensity/cgmarker/"
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
library(methyutils)
load(O)
chr <- unique(cg.mtbr$chrom)
library(BSgenome.Hsapiens.UCSC.hg38)
dna.seq <- Hsapiens[[chr]]
cposition <- GetCcontextPosition(dna.seq)
a<-cposition
b<-data.frame(1:length(a),a)
b$a[b$a==TRUE]<-1
c<-b[which(a!=0),]
d<-c
d[,2]<- -d[,2]
name<-paste("variableStep chrom=",z," span=1",sep="")
Q<-paste(file_path2,"cgmarker_n.wig",sep="")
W<-paste(file_path2,"cgmarker_p.wig",sep="")
write.table(name,file=Q,append=T,quote=F,row.names=F,col.names=F)
write.table(name,file=W,append=T,quote=F,row.names=F,col.names=F)
write.table(d,file=Q,append=T,quote=F,row.names=F,col.names=F)
write.table(c,file=W,append=T,quote=F,row.names=F,col.names=F)
}






