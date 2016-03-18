dataname<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
file_path1="/home/qzzh/breastcancer/HCC1954/HCC1954.mtbr.cg/"
#file_path2="/home/qzzh/cgDensity/coverage/RV_01/"
output=data.frame()
for (z in dataname){
O <- paste(file_path1,z,".Rdata",sep="")
load(O)
coverage<-cg.mtbr$rC_p + cg.mtbr$rC_n + cg.mtbr$rT_p + cg.mtbr$rT_n
methylevel<-with( cg.mtbr , (rC_n+rC_p)/(rC_n+rC_p+rT_n+rT_p))
methylevel[is.na(methylevel)] <- 0
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
tmp<-data.frame(chrom,posi,coverage,methylevel)
output<-rbind(output,tmp)
#Q<-paste(file_path2,"coverage.Rdata",sep="")
#write.table(outdata,file=Q,append=T,quote=F,row.names=F,col.names=F)
}
#save(output,file=Q)
library(ggplot2)
pdf1<-ggplot(output,aes(x=coverage))+geom_histogram(binwidth=1)+scale_x_continuous(limits=c(0,100),breaks=seq(0,100,5))+labs(title="HCC1954_coverage")
ggsave(pdf1,file="/home/qzzh/breastcancer/coverage_methylevel_plot/HCC1954_coverage.pdf")
pdf2<-ggplot(output,aes(x=methylevel))+geom_histogram(binwidth=0.01)+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.05))+labs(title="HCC1954_methylevel")+scale_y_continuous(limits=c(0,3000000),breaks=seq(0,3000000,500000))
ggsave(pdf2,file="/home/qzzh/breastcancer/coverage_methylevel_plot/HCC1954_methylevel.pdf")

