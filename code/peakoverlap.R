library(IRanges)
library(methyutils)
library(BSgenome.Hsapiens.UCSC.hg38)
library("data.table")
HMEC_H3k36me3 <- read.table("HMEC_H3k36me3_peaks.xls",header=T,sep="\t")
HMEC_H3k27me3 <- read.table("HMEC_H3k27me3_peaks.xls",header=T,sep="\t")
HMEC_H3k27ac <- read.table("HMEC_H3k27ac_peaks.xls",header=T,sep="\t")
HMEC_H3k9me3 <- read.table("HMEC_H3k9me3_peaks.xls",header=T,sep="\t")
HMEC_H3k4me3 <- read.table("HMEC_H3k4me3_peaks.xls",header=T,sep="\t")
HMEC_H3k4me1 <- read.table("HMEC_H3k4me1_peaks.xls",header=T,sep="\t")
HCC_H3k36me3 <- read.table("HCC_H3k36me3_peaks.xls",header=T,sep="\t")
HCC_H3k27me3 <- read.table("HCC_H3k27me3_peaks.xls",header=T,sep="\t")
HCC_H3k27ac <- read.table("HCC_H3k27ac_peaks.xls",header=T,sep="\t")
HCC_H3k9me3 <- read.table("HCC_H3k9me3_peaks.xls",header=T,sep="\t")
HCC_H3k4me3 <- read.table("HCC_H3k4me3_peaks.xls",header=T,sep="\t")
HCC_H3k4me1 <- read.table("HCC_H3k4me1_peaks.xls",header=T,sep="\t")
dmrregion<- read.table("/home/qzzh/Renbing/DMR/tmp/allregion_unique.xls",header=T,sep=",")
dna.seq <-Hsapiens[["chr1"]]
allmark.chr1 <- logical(length(dna.seq))
hmec.h3k36me3<-HMEC_H3k36me3[,1:3]
hmec.h3k27me3<-HMEC_H3k27me3[,1:3]
hmec.h3k27ac<-HMEC_H3k27ac[,1:3]
hmec.h3k9me3<-HMEC_H3k9me3[,1:3]
hmec.h3k4me3<-HMEC_H3k4me3[,1:3]
hmec.h3k4me1<-HMEC_H3k4me1[,1:3]
hcc.h3k36me3<-HCC_H3k36me3[,1:3]
hcc.h3k27me3<-HCC_H3k27me3[,1:3]
hcc.h3k27ac<-HCC_H3k27ac[,1:3]
hcc.h3k9me3<-HCC_H3k9me3[,1:3]
hcc.h3k4me3<-HCC_H3k4me3[,1:3]
hcc.h3k4me1<-HCC_H3k4me1[,1:3]
dmrregion.chr1<-dmrregion[dmrregion$chrom == "chr1",]
dmrregion.hyper.chr1<-dmrregion.chr1[dmrregion.chr1$methchangedirection == "Hyper",]
dmrregion.hypo.chr1<-dmrregion.chr1[dmrregion.chr1$methchangedirection == "Hypo",]
hmec.h3k36me3.chr1<-hmec.h3k36me3[hmec.h3k36me3$chr == "chr1",]
hmec.h3k27me3.chr1<-hmec.h3k27me3[hmec.h3k27me3$chr == "chr1",]
hmec.h3k27ac.chr1<-hmec.h3k27ac[hmec.h3k27ac$chr == "chr1",]
hmec.h3k9me3.chr1<-hmec.h3k9me3[hmec.h3k9me3$chr == "chr1",]
hmec.h3k4me3.chr1<-hmec.h3k4me3[hmec.h3k4me3$chr == "chr1",]
hmec.h3k4me1.chr1<-hmec.h3k4me1[hmec.h3k4me1$chr == "chr1",]
hcc.h3k36me3.chr1<-hcc.h3k36me3[hcc.h3k36me3$chr == "chr1",]
hcc.h3k27me3.chr1<-hcc.h3k27me3[hcc.h3k27me3$chr == "chr1",]
hcc.h3k27ac.chr1<-hcc.h3k27ac[hcc.h3k27ac$chr == "chr1",]
hcc.h3k9me3.chr1<-hcc.h3k9me3[hcc.h3k9me3$chr == "chr1",]
hcc.h3k4me3.chr1<-hcc.h3k4me3[hcc.h3k4me3$chr == "chr1",]
hcc.h3k4me1.chr1<-hcc.h3k4me1[hcc.h3k4me1$chr == "chr1",]
hmec.h3k36me3.chr1$marker <- c(1:nrow(hmec.h3k36me3.chr1))
hmec.h3k27me3.chr1$marker <- c(1:nrow(hmec.h3k27me3.chr1))
hmec.h3k27ac.chr1$marker <- c(1:nrow(hmec.h3k27ac.chr1))
hmec.h3k9me3.chr1$marker <- c(1:nrow(hmec.h3k9me3.chr1))
hmec.h3k4me3.chr1$marker <- c(1:nrow(hmec.h3k4me3.chr1))
hmec.h3k4me1.chr1$marker <- c(1:nrow(hmec.h3k4me1.chr1))
hcc.h3k36me3.chr1$marker <- c(1:nrow(hcc.h3k36me3.chr1))
hcc.h3k27me3.chr1$marker <- c(1:nrow(hcc.h3k27me3.chr1))
hcc.h3k27ac.chr1$marker <- c(1:nrow(hcc.h3k27ac.chr1))
hcc.h3k9me3.chr1$marker <- c(1:nrow(hcc.h3k9me3.chr1))
hcc.h3k4me3.chr1$marker <- c(1:nrow(hcc.h3k4me3.chr1))
hcc.h3k4me1.chr1$marker <- c(1:nrow(hcc.h3k4me1.chr1))
hmec.h3k36me3.chr1<- as.data.table(hmec.h3k36me3.chr1)
hmec.h3k27me3.chr1<- as.data.table(hmec.h3k27me3.chr1)
hmec.h3k27ac.chr1<- as.data.table(hmec.h3k27ac.chr1)
hmec.h3k9me3.chr1<- as.data.table(hmec.h3k9me3.chr1)
hmec.h3k4me3.chr1<- as.data.table(hmec.h3k4me3.chr1)
hmec.h3k4me1.chr1<- as.data.table(hmec.h3k4me1.chr1)
hcc.h3k36me3.chr1<- as.data.table(hcc.h3k36me3.chr1)
hcc.h3k27me3.chr1<- as.data.table(hcc.h3k27me3.chr1)
hcc.h3k27ac.chr1<- as.data.table(hcc.h3k27ac.chr1)
hcc.h3k9me3.chr1<- as.data.table(hcc.h3k9me3.chr1)
hcc.h3k4me3.chr1<- as.data.table(hcc.h3k4me3.chr1)
hcc.h3k4me1.chr1<- as.data.table(hcc.h3k4me1.chr1)
hmec.h3k36me3.chr1.split <- hmec.h3k36me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hmec.h3k27me3.chr1.split <- hmec.h3k27me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hmec.h3k27ac.chr1.split <- hmec.h3k27ac.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hmec.h3k9me3.chr1.split <- hmec.h3k9me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hmec.h3k4me3.chr1.split <- hmec.h3k4me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hmec.h3k4me1.chr1.split <- hmec.h3k4me1.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k36me3.chr1.split <- hcc.h3k36me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k27me3.chr1.split <- hcc.h3k27me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k27ac.chr1.split <- hcc.h3k27ac.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k9me3.chr1.split <- hcc.h3k9me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k4me3.chr1.split <- hcc.h3k4me3.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
hcc.h3k4me1.chr1.split <- hcc.h3k4me1.chr1[,.(chrom = chr, posi = c(start : end)),by = marker]
allmark.chr1[hmec.h3k36me3.chr1.split$posi]<-TRUE
allmark.chr1[hmec.h3k27me3.chr1.split$posi]<-TRUE
allmark.chr1[hmec.h3k27ac.chr1.split$posi]<-TRUE
allmark.chr1[hmec.h3k9me3.chr1.split$posi]<-TRUE
allmark.chr1[hmec.h3k4me3.chr1.split$posi]<-TRUE
allmark.chr1[hmec.h3k4me1.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k36me3.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k27me3.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k27ac.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k9me3.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k4me3.chr1.split$posi]<-TRUE
allmark.chr1[hcc.h3k4me1.chr1.split$posi]<-TRUE
rle.chr1.region <- regionAsBed(allmark.chr1,0,NULL,"chr1")
rle.chr1.IR<- IRanges(start=rle.chr1.region$start,end=rle.chr1.region$end)
hyper.chr1.IR <- IRanges(start=dmrregion.hyper.chr1$start,end=dmrregion.hyper.chr1$end)
hypo.chr1.IR <- IRanges(start=dmrregion.hypo.chr1$start,end=dmrregion.hypo.chr1$end)
hmec.h3k36me3.chr1.IR<-IRanges(start=hmec.h3k36me3.chr1$start,end=hmec.h3k36me3.chr1$end)
hmec.h3k27me3.chr1.IR<-IRanges(start=hmec.h3k27me3.chr1$start,end=hmec.h3k27me3.chr1$end)
hmec.h3k27ac.chr1.IR<-IRanges(start=hmec.h3k27ac.chr1$start,end=hmec.h3k27ac.chr1$end)
hmec.h3k9me3.chr1.IR<-IRanges(start=hmec.h3k9me3.chr1$start,end=hmec.h3k9me3.chr1$end)
hmec.h3k4me3.chr1.IR<-IRanges(start=hmec.h3k4me3.chr1$start,end=hmec.h3k4me3.chr1$end)
hmec.h3k4me1.chr1.IR<-IRanges(start=hmec.h3k4me1.chr1$start,end=hmec.h3k4me1.chr1$end)
hcc.h3k36me3.chr1.IR<-IRanges(start=hcc.h3k36me3.chr1$start,end=hcc.h3k36me3.chr1$end)
hcc.h3k27me3.chr1.IR<-IRanges(start=hcc.h3k27me3.chr1$start,end=hcc.h3k27me3.chr1$end)
hcc.h3k27ac.chr1.IR<-IRanges(start=hcc.h3k27ac.chr1$start,end=hcc.h3k27ac.chr1$end)
hcc.h3k9me3.chr1.IR<-IRanges(start=hcc.h3k9me3.chr1$start,end=hcc.h3k9me3.chr1$end)
hcc.h3k4me3.chr1.IR<-IRanges(start=hcc.h3k4me3.chr1$start,end=hcc.h3k4me3.chr1$end)
hcc.h3k4me1.chr1.IR<-IRanges(start=hcc.h3k4me1.chr1$start,end=hcc.h3k4me1.chr1$end)
ol.hyper.chr1 <- findOverlaps(rle.chr1.IR,hyper.chr1.IR)
ol.hypo.chr1 <- findOverlaps(rle.chr1.IR,hypo.chr1.IR)
ol.hmec.h3k36me3.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k36me3.chr1.IR)
ol.hmec.h3k27me3.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k27me3.chr1.IR)
ol.hmec.h3k27ac.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k27ac.chr1.IR)
ol.hmec.h3k9me3.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k9me3.chr1.IR)
ol.hmec.h3k4me3.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k4me3.chr1.IR)
ol.hmec.h3k4me1.chr1 <- findOverlaps(rle.chr1.IR,hmec.h3k4me1.chr1.IR)
ol.hcc.h3k36me3.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k36me3.chr1.IR)
ol.hcc.h3k27me3.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k27me3.chr1.IR)
ol.hcc.h3k27ac.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k27ac.chr1.IR)
ol.hcc.h3k9me3.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k9me3.chr1.IR)
ol.hcc.h3k4me3.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k4me3.chr1.IR)
ol.hcc.h3k4me1.chr1 <- findOverlaps(rle.chr1.IR,hcc.h3k4me1.chr1.IR)
chr1 <- data.frame(chr=rep("chr1",nrow(rle.chr1.region)),start=rle.chr1.region$start,end=rle.chr1.region$end,width=rle.chr1.IR@width)
chr1$HMEC_H3k36me3<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k27ac<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k27me3<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k9me3<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k4me3<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k4me1<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k36me3<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k27ac<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k27me3<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k9me3<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k4me3<-rep(0,nrow(rle.chr1.region))
chr1$HCC_H3k4me1<-rep(0,nrow(rle.chr1.region))
chr1$DMR_Hyper<-rep(0,nrow(rle.chr1.region))
chr1$DMR_Hypo<-rep(0,nrow(rle.chr1.region))
chr1$HMEC_H3k36me3[ol.hmec.h3k36me3.chr1@queryHits]<-1
chr1$HMEC_H3k27ac[ol.hmec.h3k27ac.chr1@queryHits]<-1
chr1$HMEC_H3k27me3[ol.hmec.h3k27me3.chr1@queryHits]<-1
chr1$HMEC_H3k9me3[ol.hmec.h3k9me3.chr1@queryHits]<-1
chr1$HMEC_H3k4me3[ol.hmec.h3k4me3.chr1@queryHits]<-1
chr1$HMEC_H3k4me1[ol.hmec.h3k4me1.chr1@queryHits]<-1
chr1$HCC_H3k36me3[ol.hcc.h3k36me3.chr1@queryHits]<-1
chr1$HCC_H3k27ac[ol.hcc.h3k27ac.chr1@queryHits]<-1
chr1$HCC_H3k27me3[ol.hcc.h3k27me3.chr1@queryHits]<-1
chr1$HCC_H3k9me3[ol.hcc.h3k9me3.chr1@queryHits]<-1
chr1$HCC_H3k4me3[ol.hcc.h3k4me3.chr1@queryHits]<-1
chr1$HCC_H3k4me1[ol.hcc.h3k4me1.chr1@queryHits]<-1
chr1$DMR_Hyper[ol.hyper.chr1@queryHits]<-1
chr1$DMR_Hypo[ol.hypo.chr1@queryHits]<-1
load("/home/qzzh/Renbing/DNA_methylC/HMEC/HMEC.mtbr.cg/chr1.Rdata") 
cposi <- GetCcontextPosition(dna.seq)
cg.mtbr$rC <- cg.mtbr$rC_n + cg.mtbr$rC_p
cg.mtbr$rT <- cg.mtbr$rT_n + cg.mtbr$rT_p
rC <- integer(length(dna.seq))
rC[cg.mtbr$posi] <- cg.mtbr$rC
rT <- integer(length(dna.seq))
rT[cg.mtbr$posi] <- cg.mtbr$rT
chr1$cgNum <- apply(chr1[,c("start","end")],1,function(chr1) sum(cposi[chr1["start"]:chr1["end"]]))
chr1$HMEC_rC <- apply(chr1[, c("start", "end")], 1, function(chr1) sum(rC[chr1["start"]:chr1["end"]]))
chr1$HMEC_rT <- apply(chr1[, c("start", "end")], 1, function(chr1) sum(rT[chr1["start"]:chr1["end"]]))
chr1$HMEC_score <- chr1$HMEC_rC/(chr1$HMEC_rC + chr1$HMEC_rT)
load("/home/qzzh/Renbing/DNA_methylC/HCC1954/HCC1954.mtbr.cg/chr1.Rdata")
cg.mtbr$rC <- cg.mtbr$rC_n + cg.mtbr$rC_p
cg.mtbr$rT <- cg.mtbr$rT_n + cg.mtbr$rT_p
rC[cg.mtbr$posi] <- cg.mtbr$rC
rT[cg.mtbr$posi] <- cg.mtbr$rT
chr1$HCC_rC <- apply(chr1[, c("start", "end")], 1, function(chr1) sum(rC[chr1["start"]:chr1["end"]]))
chr1$HCC_rT <- apply(chr1[, c("start", "end")], 1, function(chr1) sum(rT[chr1["start"]:chr1["end"]]))
chr1$HCC_score <- chr1$HCC_rC/(chr1$HCC_rC + chr1$HCC_rT)
save(chr1,file="/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/chr1.Rdata")
