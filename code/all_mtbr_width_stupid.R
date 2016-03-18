
file_path2<-"/home/qzzh/J.R.Ecker/all_mtbr_width/"
load("/home/qzzh/J.R.Ecker/mtbr/AD_02/AD_02_unique.mtbr.cg/chrY.Rdata")
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
AD_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
AD_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
AD_02_score<-AD_02_rC/(AD_02_rC+AD_02_rT)
AD_02_score[is.na(AD_02_score[])] <- 0
AD_02<-data.frame(chrom,posi,AD_02_rC,AD_02_rT,AD_02_score)
load("/home/qzzh/J.R.Ecker/mtbr/EG_02/EG_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
EG_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
EG_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
EG_02_score<-EG_02_rC/(EG_02_rC+EG_02_rT)
EG_02_score[is.na(EG_02_score[])] <- 0
EG_02<-data.frame(posi,EG_02_rC,EG_02_rT,EG_02_score)
all<-merge(AD_02,EG_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/FT_02/FT_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
FT_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
FT_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
FT_02_score<-FT_02_rC/(FT_02_rC+FT_02_rT)
FT_02_score[is.na(FT_02_score[])] <- 0
FT_02<-data.frame(posi,FT_02_rC,FT_02_rT,FT_02_score)
all<-merge(all,FT_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/LG_01/LG_01_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
LG_01_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
LG_01_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
LG_01_score<-LG_01_rC/(LG_01_rC+LG_01_rT)
LG_01_score[is.na(LG_01_score[])] <- 0
LG_01<-data.frame(posi,LG_01_rC,LG_01_rT,LG_01_score)
all<-merge(all,LG_01,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/PA_02/PA_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
PA_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
PA_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
PA_02_score<-PA_02_rC/(PA_02_rC+PA_02_rT)
PA_02_score[is.na(PA_02_score[])] <- 0
PA_02<-data.frame(posi,PA_02_rC,PA_02_rT,PA_02_score)
all<-merge(all,PA_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/PO_02/PO_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
PO_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
PO_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
PO_02_score<-PO_02_rC/(PO_02_rC+PO_02_rT)
PO_02_score[is.na(PO_02_score[])] <- 0
PO_02<-data.frame(posi,PO_02_rC,PO_02_rT,PO_02_score)
all<-merge(all,PO_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/RV_01/RV_01_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
RV_01_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
RV_01_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
RV_01_score<-RV_01_rC/(RV_01_rC+RV_01_rT)
RV_01_score[is.na(RV_01_score[])] <- 0
RV_01<-data.frame(posi,RV_01_rC,RV_01_rT,RV_01_score)
all<-merge(all,RV_01,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/SB_02/SB_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
SB_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
SB_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
SB_02_score<-SB_02_rC/(SB_02_rC+SB_02_rT)
SB_02_score[is.na(SB_02_score[])] <- 0
SB_02<-data.frame(posi,SB_02_rC,SB_02_rT,SB_02_score)
all<-merge(all,SB_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/SB_03/SB_03_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
SB_03_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
SB_03_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
SB_03_score<-SB_03_rC/(SB_03_rC+SB_03_rT)
SB_03_score[is.na(SB_03_score[])] <- 0
SB_03<-data.frame(posi,SB_03_rC,SB_03_rT,SB_03_score)
all<-merge(all,SB_03,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/AO_02/AO_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
AO_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
AO_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
AO_02_score<-AO_02_rC/(AO_02_rC+AO_02_rT)
AO_02_score[is.na(AO_02_score[])] <- 0
AO_02<-data.frame(posi,AO_02_rC,AO_02_rT,AO_02_score)
all<-merge(all,AO_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/FT_01/FT_01_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
FT_01_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
FT_01_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
FT_01_score<-FT_01_rC/(FT_01_rC+FT_01_rT)
FT_01_score[is.na(FT_01_score[])] <- 0
FT_01<-data.frame(posi,FT_01_rC,FT_01_rT,FT_01_score)
all<-merge(all,FT_01,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/OV_02/OV_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
OV_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
OV_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
OV_02_score<-OV_02_rC/(OV_02_rC+OV_02_rT)
OV_02_score[is.na(OV_02_score[])] <- 0
OV_02<-data.frame(posi,OV_02_rC,OV_02_rT,OV_02_score)
all<-merge(all,OV_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/GA_02/GA_02_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
GA_02_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
GA_02_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
GA_02_score<-GA_02_rC/(GA_02_rC+GA_02_rT)
GA_02_score[is.na(GA_02_score[])] <- 0
GA_02<-data.frame(posi,GA_02_rC,GA_02_rT,GA_02_score)
all<-merge(all,GA_02,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/PO_01/PO_01_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
PO_01_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
PO_01_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
PO_01_score<-PO_01_rC/(PO_01_rC+PO_01_rT)
PO_01_score[is.na(PO_01_score[])] <- 0
PO_01<-data.frame(posi,PO_01_rC,PO_01_rT,PO_01_score)
all<-merge(all,PO_01,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/SX_01/SX_01_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
SX_01_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
SX_01_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
SX_01_score<-SX_01_rC/(SX_01_rC+SX_01_rT)
SX_01_score[is.na(SX_01_score[])] <- 0
SX_01<-data.frame(posi,SX_01_rC,SX_01_rT,SX_01_score)
all<-merge(all,SX_01,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)
load("/home/qzzh/J.R.Ecker/mtbr/RA_03/RA_03_unique.mtbr.cg/chrY.Rdata")
posi<-cg.mtbr$posi
RA_03_rC<-cg.mtbr$rC_n + cg.mtbr$rC_p
RA_03_rT<-cg.mtbr$rT_n + cg.mtbr$rT_p
RA_03_score<-RA_03_rC/(RA_03_rC+RA_03_rT)
RA_03_score[is.na(RA_03_score[])] <- 0
RA_03<-data.frame(posi,RA_03_rC,RA_03_rT,RA_03_score)
all<-merge(all,RA_03,by.x="posi",by.y="posi",all.x=TRUE,all.y=TRUE)

all$chrom<-"chrY"
all[is.na(all[])]<-0
#library("matrixStats")
#s<-grep("score",colnames(all))
#score<-all[,s]
#score[,1]<-as.numeric(score[,1])
#score[,2]<-as.numeric(score[,2])
#score[,3]<-as.numeric(score[,3])
#score[,4]<-as.numeric(score[,4])
#score[,5]<-as.numeric(score[,5])
#score[,6]<-as.numeric(score[,6])
#score[,7]<-as.numeric(score[,7])
#score[,8]<-as.numeric(score[,8])
#score[,9]<-as.numeric(score[,9])
#score[,10]<-as.numeric(score[,10])
#score[,11]<-as.numeric(score[,11])
#score[,12]<-as.numeric(score[,12])
#score[,13]<-as.numeric(score[,13])
#score[,14]<-as.numeric(score[,14])
#score[,15]<-as.numeric(score[,15])
#score[,16]<-as.numeric(score[,16])
#scorematrix<-as.matrix(score)
#sd<-rowSds(scorematrix)
#all$sd<-sd
O<-paste(file_path2,"all_mtbr_width_chrY.Rdata",sep="")
save(all,file=O)

