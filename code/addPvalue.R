library("matrixStats")
library("data.table")
library("reshape2")
library("methyutils")
MultiAddPv <- function(mg.data){
                   sn <- (ncol(mg.data)-2) / 2
                   mg.data <- mg.data[order(mg.data$chrom,mg.data$posi),]
                   m <- melt(mg.data, id=c("chrom", "posi"),variable.name="read", value.name="count")
                   m <- m[order(m$chrom,m$posi),]
                   cp <- MultiChisqTest(m$count,nrow=sn,ncol=2)
                   mg.data$pv<-cp
                   return(mg.data)
                                }
file_path1 <- "/home/qzzh/J.R.Ecker/all_mtbr_width/"
file_path2 <- "/home/qzzh/J.R.Ecker/addPvalue/"
data <- list.files("/home/qzzh/J.R.Ecker/all_mtbr_width/")
for (i in data){
I <- paste(file_path1,i,sep="")
load(I)
r<- grep("rC|rT",colnames(output))
read <- output[,r]
coverage <- rowSums(read)
chrom <- output$chrom
posi <- output$posi
test <- cbind(chrom,posi,read)   
addPvalue <- MultiAddPv(test)
s <- grep("score",colnames(output))
score <- output[,s]
mean <- rowMeans(score)
scorematrix <- as.matrix(score)
sd <- rowSds(scorematrix)
addPvalue$Score_sd <- sd
addPvalue$coverage <- coverage
addPvalue$Score_mean <- mean
addPvalue[addPvalue < 0] <- 1
O <- paste(file_path2,i,sep="")
save(addPvalue,file=O)
               }

