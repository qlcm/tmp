file<-list.files("/home/qzzh/J.R.Ecker/all_mtbr_width/")
file_path<-"/home/qzzh/J.R.Ecker/all_mtbr_width/"
output<-"/home/qzzh/J.R.Ecker/change_all_mtbr_width/"
for(i in file){
I<-paste(file_path,i,sep="")
load(I)
chrom<-all$chrom
posi<-all$posi
r <- grep("rC|rT", colnames(all)) 
read <- all[,r]
all<-cbind(chrom,posi,read)
O<-paste(output,i,sep="")
save(all,file=O)
}


