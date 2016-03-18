dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
file_path1="/home/qzzh/J.R.Ecker/mtbr/SB_02/SB_02_unique.mtbr.cg/"
file_path2="/home/qzzh/cgDensity/gap/SB_02/"
file_path3="/home/qzzh/cgDensity/overlap/SB_02/"
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
library(methyutils)
load(O)
output<-cgDensity(cg.mtbr,genome = "hg38",window = 1250)
P<-paste(file_path3,z,sep="")
Q<-paste(file_path2,z,sep="")
write.table(output$gap,file=Q,append=F,quote=F,row.names=F,col.names=F)
write.table(output$overlap,file=P,append=F,quote=F,row.names=F,col.names=F)
message(z," is running finished ",date())
}






