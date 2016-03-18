library(ggplot2)
load("/home/qzzh/breastcancer/HMEC/HMEC.mtbr.cg/chr19.Rdata")
HMECmeth <- with( cg.mtbr , (rC_n+rC_p)/(rC_n+rC_p+rT_n+rT_p))
HMECmeth[is.na(HMECmeth)] <- 0
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
HMEC <- data.frame(chrom,posi,HMECmeth)
load("/home/qzzh/breastcancer/HCC1954/HCC1954.mtbr.cg/chr19.Rdata")
HCCmeth <- with( cg.mtbr , (rC_n+rC_p)/(rC_n+rC_p+rT_n+rT_p))
HCCmeth[is.na(HCCmeth)] <- 0
chrom<-cg.mtbr$chrom
posi<-cg.mtbr$posi
HCC <- data.frame(chrom,posi,HCCmeth)
all <- merge(HMEC,HCC,all=T)
all[is.na(all)] <- 0
all$chrom <- "chr19"
save(all,file="/home/qzzh/breastcancer/Rdata/HMEC-HCC.Rdata")
p <- ggplot(all,aes(x=all$HCCmeth,y=all$HMECmeth))+geom_point(color="blue",alpha=0.1)+labs(x="HCC1954",y="HMEC")
ggsave(p,file="/home/qzzh/breastcancer/coverage_methylevel_plot/HMEC-HCC.pdf")
