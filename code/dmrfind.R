library(methyutils)
file_path1 <- "/home/qzzh/J.R.Ecker/addPvalue/"
file_path2 <- "/home/qzzh/J.R.Ecker/dmrFinder/"
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
           "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
dmrFinder <- function(mg.data, pvbottom=0, pvtop=0.05, cf.length=10, tolerance.length=2){
                     chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                              "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
                      dmr <- data.frame(chrom = as.character(), start = as.integer(), end = as.integer())
                               for(chr in chrs){
                                   chr.mg <- mg.data[mg.data$chrom == chr,]
                                   dl <- chr.mg$pv > pvbottom & chr.mg$pv <= pvtop
                                   region.index  <- regionAsBed(dl, cf.length, tolerance.length=2,chrom=chr)
                                   region.postion <- data.frame(chrom   = as.character(chr.mg$chrom[region.index$start]),
                                   start = as.integer(chr.mg$posi[region.index$start]),
                                   end   = as.integer(chr.mg$posi[region.index$end]))
                                   dmr <- data.frame(chrom = c(as.character(dmr$chrom),as.character(region.postion$chrom)),
                                   start = c(dmr$start,region.postion$start),
                                   end   = c(dmr$end,region.postion$end))
                                          }

                           return(dmr)

                                    }
for (i in chrs) {
          I <- paste(file_path1,"all_mtbr_width_",i,".Rdata",sep="")
          load(I)
          dmrregion <- dmrFinder(addPvalue,pvbottom=0, pvtop=0.05, cf.length=10, tolerance.length=2)
          O <- paste(file_path2,i,sep="")
          write.table(dmrregion,file=O,append=F,quote=F,row.names=F,col.names=F)
          message(i," is running finished ",date())
            }
     

