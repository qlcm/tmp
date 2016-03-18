library("data.table")
load("/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/chrX.Rdata")
df  <- df[order(df$start),]
data <- c("HMEC_H3k36me3","HCC_H3k36me3","HMEC_H3k27ac","HCC_H3k27ac","HMEC_H3k27me3","HCC_H3k27me3")
for (i in data){
I <- paste("/home/qzzh/Renbing/Chip-seq/ToGB/",i,".bdg",sep="")
tmp1 <- fread(i)
tmp2 <- tmp1[tmp1$V1 == "chrX",]
colnames(tmp2) <- c("chrom","start","end","value")
tmp2 <- tmp2[tmp2$value > 0 ,]
tmp2$start <- tmp2$start + 1
tmp2$marker <- c(1:nrow(tmp2))
dt.tmp2 <- tmp2[,.(chrom = chrom, posi = c(start : end), value = value),by = marker]
dt.tmp2$marker <- NULL
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
tmp3 <- data.frame(value=df.add.value$score)
colnames(tmp3) <- c(i)
df <- cbind(df,tmp3)
}
tmp <- grep("[^df]",ls())
rm(list=ls()[tmp])
gc()

save(df,file="/home/qzzh/Renbing/Chip-seq/macsresult/peakRdata/chrX.Rdata")
