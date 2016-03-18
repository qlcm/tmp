tissuename<-c("AD_02","EG_02","FT_02","LG_01","PA_02","PO_02","RV_01","SB_03","AO_02","FT_01","GA_02","OV_02","PO_01","RA_03","SB_02","SX_01")
file_path1<-"/home/qzzh/J.R.Ecker/mtbr/"
file_path2<-"/home/qzzh/J.R.Ecker/all_mtbr/"
output<-data.frame()
for (z in tissuename){
I<-paste(file_path1,z,"/",z,"_unique.mtbr.cg","/","chrY.Rdata",sep="")
load(I)
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
cg.mtbr$tissue<-z
tissue<-cg.mtbr$tissue
tmp<-data.frame(chrom,posi,rC,rT,tissue)
output<-rbind(output,tmp)
}
O<-paste(file_path2,"all_mtbr_long_chrY.Rdata",sep="")
save(output,file=O)

