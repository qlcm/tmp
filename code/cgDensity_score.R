dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
file_path1="/home/qzzh/J.R.Ecker/PO_02_unique.mtbr.cg/"
file_path2="/home/qzzh/J.R.Ecker/cgDensity/use_hg38/score/PO_02/"
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
library(methyutils)
load(O)
cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
chr <- unique(cg.mtbr$chrom)
library(BSgenome.Hsapiens.UCSC.hg38)
dna.seq <- Hsapiens[[chr]]
rC <- integer(length(dna.seq))
rC[cg.mtbr$posi] <- cg.mtbr$rC
rT <- integer(length(dna.seq))
rT[cg.mtbr$posi] <- cg.mtbr$rT
win <- list(L = 500, R = 500)
name<-paste("variableStep chrom=",z," span=1",sep="")
rCs <- swsCalc(rC, win)
rTs <- swsCalc(rT, win)
score <- rCs/(rCs + rTs)
score[is.na(score[])] <- 0
SCORE<-data.frame(1:length(score),score)
score_p<-SCORE[which(score!=0),]
score_n<-score_p
score_n[,2]<- -score_n[,2]
S<-paste(file_path2,"score_n.wig",sep="")
P<-paste(file_path2,"score_p.wig",sep="")
write.table(name,file=S,append=T,quote=F,row.names=F,col.names=F)
write.table(name,file=P,append=T,quote=F,row.names=F,col.names=F)
write.table(score_n,file=S,append=T,quote=F,row.names=F,col.names=F)
write.table(score_p,file=P,append=T,quote=F,row.names=F,col.names=F)
}






