dataname<-c("chr22")
filenames<-list.files("/home/qzzh/J.R.Ecker/mtbr")
file_path1="/home/qzzh/J.R.Ecker/mtbr/"
file_path2="/home/qzzh/J.R.Ecker/all_mtbr_width/"
output=data.frame()
load("/home/qzzh/J.R.Ecker/mtbr/AD_02/AD_02_unique.mtbr.cg/chrY.Rdata")
chrom <- cg.mtbr$chrom
posi<-cg.mtbr$posi
rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
score <- rC / (rC+rT)
score[is.na(score[])] <- 0
test<-data.frame(chrom,posi,rC,rT,score)
output<-merge(output,test)
output$rC<-NULL
output$rT<-NULL
for (samplename in filenames){
message(samplename," is running",date())
for (z in dataname){
O <- paste(file_path1,samplename,"/",samplename,"_unique.mtbr.cg","/",z,".Rdata",sep="")
load(O)
rC <- cg.mtbr$rC_p + cg.mtbr$rC_n 
rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
score <- rC / (rC+rT)
score[is.na(score[])] <- 0
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
tmp<-data.frame(chrom,posi,rC,rT,score)
colname1 <- paste(samplename,"_rC",sep="")
colname2 <- paste(samplename,"_rT",sep="")
colname3 <- paste(samplename,"_score",sep="")
colnames(tmp)<-c("chrom","posi",colname1,colname2,colname3)
output<-merge(output,tmp,all=T)
}
output$chrom <-z
output[is.na(output[])]<-0
P <- paste(file_path2,z,".Rdata",sep="")
save(output,file=P)
}
