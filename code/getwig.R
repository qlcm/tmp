dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
file_path1="/home/qzzh/J.R.Ecker/SB_03/SB_03_unique.mtbr.cg/"
wigfilePath="/home/qzzh/J.R.Ecker/SB_03/"
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
load(O)
head <- paste("variableStep chrom=",z, " span=1\n", sep="")
Sp <- data.frame(posi=cg.mtbr$posi,score=cg.mtbr$rC_p/(cg.mtbr$rC_p + cg.mtbr$rT_p))
Sp <- Sp[!is.na(Sp$score),]
Sp <- Sp[Sp$score>0,]
Rmp <- data.frame(posi=cg.mtbr$posi,score=cg.mtbr$rC_p)
Rmp <- Rmp[!is.na(Rmp$score),]
Rmp <- Rmp[Rmp$score>0,]
Rtotalp <- data.frame(posi=cg.mtbr$posi,score=(cg.mtbr$rC_p + cg.mtbr$rT_p))
Rtotalp <- Rtotalp[!is.na(Rtotalp$score),]
Rtotalp <- Rtotalp[Rtotalp$score>0,]
Sn <- data.frame(posi=cg.mtbr$posi+1,score=(-cg.mtbr$rC_n/(cg.mtbr$rC_n + cg.mtbr$rT_n)))
Sn <- Sn[!is.na(Sn$score),]
Sn <- Sn[Sn$score<0,]
Sn$posi <- as.integer(Sn$posi)
Rmn <- data.frame(posi=cg.mtbr$posi+1,score=-cg.mtbr$rC_n)
Rmn <- Rmn[!is.na(Rmn$score),]
Rmn <- Rmn[Rmn$score<0,]
Rmn$posi <- as.integer(Rmn$posi)
Rtotaln <- data.frame(posi=cg.mtbr$posi+1,score=-(cg.mtbr$rC_n + cg.mtbr$rT_n))
Rtotaln <- Rtotaln[!is.na(Rtotaln$score),]
Rtotaln <- Rtotaln[Rtotaln$score<0,]
Rtotaln$posi <- as.integer(Rtotaln$posi)
Spf <- paste(wigfilePath,"methy.score_p.wig",sep="")
Rmpf <- paste(wigfilePath,"methy.cover_p.wig",sep="")
Rtotalpf <- paste(wigfilePath,"total.cover_p.wig",sep="")
Rtotalnf <- paste(wigfilePath,"total.cover_n.wig",sep="")
Snf <- paste(wigfilePath,"methy.score_n.wig",sep="")
Rmnf <- paste(wigfilePath,"methy.cover_n.wig",sep="")
cat(head, file=Spf,append=T)
write.table(Sp, Spf, row.names=F, col.names=F, quote=F, append=T)
cat(head, file=Rmpf,append=T)
write.table(Rmp, Rmpf, row.names=F, col.names=F, quote=F, append=T)
cat(head, file=Rtotalpf,append=T)
write.table(Rtotalp, Rtotalpf, row.names=F, col.names=F, quote=F, append=T)
cat(head, file=Snf,append=T)
write.table(Sn, Snf, row.names=F, col.names=F, quote=F, append=T)
cat(head, file=Rmnf,append=T)
write.table(Rmn, Rmnf, row.names=F, col.names=F, quote=F, append=T)
cat(head, file=Rtotalnf,append=T)
write.table(Rtotaln, Rtotalnf, row.names=F, col.names=F, quote=F, append=T)
}










