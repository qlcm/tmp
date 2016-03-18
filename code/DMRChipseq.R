library("data.table")
load("/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/chrX.Rdata")
chrX <- chrX[order(chrX$start),]
HMEC_H3k4me1 <- fread("/home/qzzh/Renbing/Chip-seq/ToGB/HMEC_H3k4me1.bdg")
HCC_H3k4me1 <- fread("/home/qzzh/Renbing/Chip-seq/ToGB/HCC_H3k4me1.bdg")
chrX.HMEC.H3k4me1 <- HMEC_H3k4me1[HMEC_H3k4me1$V1 == "chrX",]
chrX.HCC.H3k4me1 <- HCC_H3k4me1[HCC_H3k4me1$V1 == "chrX",]
colnames(chrX.HMEC.H3k4me1 ) <- c("chrom","start","end","value")
colnames(chrX.HCC.H3k4me1 ) <- c("chrom","start","end","value")
chrX.HMEC.H3k4me1 <- chrX.HMEC.H3k4me1[chrX.HMEC.H3k4me1$value > 0 ,]
chrX.HCC.H3k4me1 <- chrX.HCC.H3k4me1[chrX.HCC.H3k4me1$value > 0 ,]
chrX.HMEC.H3k4me1$start <- chrX.HMEC.H3k4me1$start + 1
chrX.HCC.H3k4me1$start <- chrX.HCC.H3k4me1$start + 1
chrX.HMEC.H3k4me1$marker <- c(1:nrow(chrX.HMEC.H3k4me1))
chrX.HCC.H3k4me1$marker <- c(1:nrow(chrX.HCC.H3k4me1))
dt.chrX.HMEC.H3k4me1 <- chrX.HMEC.H3k4me1[,.(chrom = chrom, posi = c(start : end), value = value),by = marker]
dt.chrX.HCC.H3k4me1 <- chrX.HCC.H3k4me1[,.(chrom = chrom, posi = c(start : end), value = value),by = marker]
dt.chrX.HMEC.H3k4me1$marker <- NULL
dt.chrX.HCC.H3k4me1$marker <- NULL
chrX.region <- chrX[,1:3]
colnames(chrX.region) <-  c("chrom","start","end")
chrX.region$marker <- c(1:nrow(chrX.region))
dt.chrX.region <- as.data.table(chrX.region)
dt.chrX.region$chrom <- as.character(dt.chrX.region$chrom)
dt.chrX.region$start <- as.integer(dt.chrX.region$start)
dt.chrX.region$end <- as.integer(dt.chrX.region$end)
dt.chrX.region.split <- dt.chrX.region[,.(chrom = chrom, posi = c(start : end)),by = marker]
dt.HMEC.value <- merge(dt.chrX.HMEC.H3k4me1,dt.chrX.region.split,by.x = c("posi","chrom"),by.y = c("posi","chrom"))
dt.HCC.value <- merge(dt.chrX.HCC.H3k4me1,dt.chrX.region.split,by.x = c("posi","chrom"),by.y = c("posi","chrom"))
dt.HMEC.chip.value <- dt.HMEC.value[,.(H3k4me1 = sum(value)),by = marker]
dt.HCC.chip.value <- dt.HCC.value[,.(H3k4me1 = sum(value)),by = marker]
chrX.HMEC.add.value <- merge(dt.chrX.region,dt.HMEC.chip.value,by="marker",all = TRUE)
chrX.HCC.add.value <- merge(dt.chrX.region,dt.HCC.chip.value,by="marker",all = TRUE)
chrX$HMEC_H3k4me1_value <- chrX.HMEC.add.value$H3k4me1
chrX$HCC_H3k4me1_value <- chrX.HCC.add.value$H3k4me1
tmp <- grep("[^chrX]",ls())
rm(list=ls()[tmp])
gc()

save(chrX,file="/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/chrX.Rdata")
