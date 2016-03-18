GetDensity <- function(ref.fa, kWinSize) {
	
	posi <- matchPattern("CG",ref.fa)
	rt <- logical(ref.length)
	rt[start(posi)] <- TRUE
   	win <- list(L = kWinSize, R = kWinSize)
 	return(swsCalc(rt, win))
}

GetScore <- function(cg.mtbr, kWinSize, ref.length) {

	  colnames(cg.mtbr) <- c("chr", "posi", "rC_n", "rC_p", "rT_n", "rT_p")  
	
	  cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
	  cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
	  
	  rC <- integer(ref.length)
	  rC[cg.mtbr$posi] <- cg.mtbr$rC
	  rT <- integer(ref.length)
	  rT[cg.mtbr$posi] <- cg.mtbr$rT
	  win <- list(L = kWinSize, R = kWinSize)
	  rCs <- swsCalc(rC, win)
	  rTs <- swsCalc(rT, win)
	  score <- rCs/(rCs + rTs)
	  score[is.na(score[])] <- 0
	  
  return(score)
}

library("methyutils")
library("BSgenome.Hsapiens.UCSC.hg38")
tissues <- list.files("/home/qzzh/J.R.Ecker/mtbr/mtbr2/")
file_path<- "/home/qzzh/J.R.Ecker/mtbr/mtbr2/"
#chrs <- c("chr21","chr22","chrX","chrY")
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
           "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
#ref.fa <- Hsapiens[[ts]]
#ref.length <- length(ref.fa)
#ts.list <- list()
for(ts in tissues){
for (chr in chrs){
	message(paste(ts," ",chr," ","is running ",sep=""),date())
	load(paste(file_path,ts,"/",ts,"_unique.mtbr.cg/",chr,".Rdata",sep=""))
        ref.fa <- Hsapiens[[chr]]
        ref.length <- length(ref.fa)
	#ks <- c()
	#cors <- c()		
	kWinSize <- 1250 	
		density <- GetDensity(ref.fa, kWinSize)
		score <- GetScore(cg.mtbr, kWinSize, ref.length)
		cr <- cor(density[cg.mtbr$posi],score[cg.mtbr$posi])
		#ks <- c(ks,kWinSize)
		#cors <- c(cors,cr)	
	
	sw.df <- data.frame(winsize = kWinSize,cor = cr,tissue = ts,chrom=chr)
	#ts.list[[ts]] <- sw.df
       write.table(sw.df,"/home/qzzh/cgDensity/SelectSlidingwindow/SelectSlidingwindow_mtbr2.txt",row.names=F,col.names=F,quote=F,append=T)
}
}
#tissue <- do.call(rbind,ts.list)

#write.table(tissue,"/home/qzzh/cgDensity/SelectSlidingwindow/SelectSlidingwindow.txt",row.names=F,col.names=T,quote=F)
