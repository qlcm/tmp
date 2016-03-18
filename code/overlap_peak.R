library(IRanges)
library(methyutils)
library(BSgenome.Hsapiens.UCSC.hg38)
library("data.table")
chroms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
dataname <- c("HMEC_H3k36me3","HCC_H3k36me3","HMEC_H3k27ac","HCC_H3k27ac","HMEC_H3k27me3","HCC_H3k27me3",
"HMEC_H3k9me3","HCC_H3k9me3","HMEC_H3k4me3","HCC_H3k4me3","HMEC_H3k4me1","HCC_H3k4me1") 
input1 <- "/home/qzzh/Renbing/"
input2 <- "/home/qzzh/Renbing/DNA_methylC/HMEC/HMEC.mtbr.cg/"  
input3 <- "/home/qzzh/Renbing/DNA_methylC/HCC1954/HCC1954.mtbr.cg/"  
output <- "/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/"
for (chr in chroms){
  dna.seq <-Hsapiens[[chr]]
  allmark <- logical(length(dna.seq))
    for (name in dataname){
         I <- paste(input1,name,"_peaks.xls",sep="") 
         tmp <- read.table(I,header=T,sep="\t")
         chip.peak <- tmp[,1:3]
         chip.peak.chrom <-chip.peak[chip.peak$chr == chr,]
         chip.peak.chrom$marker <- c(1:nrow(chip.value.chrom))
         chip.peak.chrom <- as.data.table(chip.value.chrom)
         chip.peak.chrom.split <- chip.peak.chrom[,.(chrom = chr, posi = c(start : end)),by = marker]
         allmark[chip.value.chrom.split$posi] <- TRUE
                          }
  rle.region <- regionAsBed(allmark,0,NULL,chr)
  rle.IR <- IRanges(start=rle.region$start,end=rle.region$end)
  df <- data.frame(chr=rep(chr,nrow(rle.region)),start=rle.region$start,end=rle.region$end,width=rle.IR@width)
    for (name in dataname){
         I <- paste(input1,name,"_peaks.xls",sep="") 
         tmp <- read.table(I,header=T,sep="\t")
         chip.peak <- tmp[,1:3]
         chip.peak.chrom <-chip.peak[chip.peak$chr == chr,]
         chip.peak.chrom.IR <-IRanges(start=chip.peak.chrom$start,end=chip.peak.chrom$end)
         ol <- findOverlaps(rle.IR,chip.value.chrom.IR)
         tmp1 <-data.frame( marker=rep(0,nrow(rle.region)))
         tmp1$marker[ol@queryHits] <- 1
         colnames(tmp1) <- c(name)
         df <- cbind(df,tmp1)
                           }
   V1 <- paste(input2,chr,".Rdata",sep="")
   load(V1) 
   cposi <- GetCcontextPosition(dna.seq)
   cg.mtbr$rC <- cg.mtbr$rC_n + cg.mtbr$rC_p
   cg.mtbr$rT <- cg.mtbr$rT_n + cg.mtbr$rT_p
   rC <- integer(length(dna.seq))
   rC[cg.mtbr$posi] <- cg.mtbr$rC
   rT <- integer(length(dna.seq))
   rT[cg.mtbr$posi] <- cg.mtbr$rT
   df$cgNum <- apply(df[,c("start","end")],1,function(df) sum(cposi[df["start"]:df["end"]]))
   df$HMEC_rC <- apply(df[, c("start", "end")], 1, function(df) sum(rC[df["start"]:df["end"]]))
   df$HMEC_rT <- apply(df[, c("start", "end")], 1, function(df) sum(rT[df["start"]:df["end"]]))
   df$HMEC_score <- df$HMEC_rC/(df$HMEC_rC + df$HMEC_rT)
   V2 <- paste(input3,chr,".Rdata",sep="")
   load(V2)
   cg.mtbr$rC <- cg.mtbr$rC_n + cg.mtbr$rC_p
   cg.mtbr$rT <- cg.mtbr$rT_n + cg.mtbr$rT_p
   rC[cg.mtbr$posi] <- cg.mtbr$rC
   rT[cg.mtbr$posi] <- cg.mtbr$rT
   df$HCC_rC <- apply(df[, c("start", "end")], 1, function(df) sum(rC[df["start"]:df["end"]]))
   df$HCC_rT <- apply(df[, c("start", "end")], 1, function(df) sum(rT[df["start"]:df["end"]]))
   df$HCC_score <- df$HCC_rC/(df$HCC_rC + df$HCC_rT)
   df  <- df[order(df$start),]
   for (name in dataname){
     I <- paste("/home/qzzh/Renbing/Chip-seq/ToGB/",name,".bdg",sep="")
     tmp2 <- fread(i)
     tmp3 <- tmp2[tmp2$V1 == chr,]
     colnames(tmp3) <- c("chrom","start","end","value")
     tmp3 <- tmp3[tmp3$value > 0 ,]
     tmp3$start <- tmp3$start + 1
     tmp3$marker <- c(1:nrow(tmp3))
     dt.tmp3 <- tmp3[,.(chrom = chrom, posi = c(start : end), value = value),by = marker]
     dt.tmp3$marker <- NULL
     df.region <- df[,1:3]
     colnames(df.region) <-  c("chrom","start","end")
     df.region$marker <- c(1:nrow(df.region))
     dt.df.region <- as.data.table(df.region)
     dt.df.region$chrom <- as.character(dt.df.region$chrom)
     dt.df.region$start <- as.integer(dt.df.region$start)
     dt.df.region$end <- as.integer(dt.df.region$end)
     dt.df.region.split <- dt.df.region[,.(chrom = chrom, posi = c(start : end)),by = marker]
     dt.value <- merge(dt.tmp2,dt.df.region.split,by.x = c("posi","chrom"),by.y = c("posi","chrom"))
     dt.chip.value <- dt.value[,.(score = sum(value)),by = marker]
     df.add.value <- merge(dt.df.region,dt.chip.value,by="marker",all = TRUE)
     tmp4 <- data.frame(value=df.add.value$score)
     colnames(tmp4) <- c(name)
     df <- cbind(df,tmp4)                
                            }
     tmp5 <- grep("[^df]",ls())
     rm(list=ls()[tmp5])
     gc()
  O <- paste(output,chr,".Rdata",seq="")
  save(df,file=O)

              }
