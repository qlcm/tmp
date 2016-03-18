library(ggplot2)
########################################################################################
################################Figure5-split###########################################
########################################################################################
plot1 <- read.table("/home/qzzh/plot/bonemarrowandKidneyExp-barplot.txt",header = T,sep="\t")
plot2 <- read.table("/home/qzzh/plot/bonemarrowandKidneyExp-lineplot.txt",header = T)
plot3 <- read.table("/home/qzzh/plot/bonemarrowandKidneyPhy-barplot.txt",header = T,sep = "\t")
plot4 <- read.table("/home/qzzh/plot/bonemarrowandKidneyPhy-lineplot.txt",header = T)
plot1 <- plot1[order(plot1$negativelog10),]
plot3 <- plot3[order(plot3$nagativelog10),]
plot1$Desc <- with(plot1, factor(Desc, levels=Desc, ordered=TRUE))
plot3$Desc <- with(plot3, factor(Desc, levels=Desc, ordered=TRUE))
p1 <- ggplot(plot1,aes(x=plot1$Desc,y=plot1$negativelog10))+theme(panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(size=0.5),legend.title=element_blank())+geom_bar(fill="blue",stat="identity")+theme(axis.text.x=element_text(angle=0,colour="black",size=15,hjust=0.5))+theme(axis.text.y=element_text(angle=0,colour="black",size=10,hjust=1))+coord_flip()+labs(x="Desc",y="negativelog10(p-value)")+scale_y_discrete(limits=c(0,50,100,150))
p2 <- ggplot(plot2,aes(x=plot2$Rank,y=plot2$percentage))+geom_point(size=3,color="blue")+ylim(c(0, 1)) + geom_line(size=1.2,color="blue")+labs(x="Rank",y="percentage")
p3 <- ggplot(plot3,aes(x=plot3$Desc,y=plot3$nagativelog10))+theme(panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(size=0.5),legend.title=element_blank())+geom_bar(fill="blue",stat="identity")+theme(axis.text.x=element_text(angle=0,colour="black",size=15,hjust=0.5))+theme(axis.text.y=element_text(angle=0,colour="black",size=10,hjust=1))+coord_flip()+labs(x="Desc",y="negativelog10(p-value)")+scale_y_discrete(limits=c(0,50,100,150,200))
p4 <- ggplot(plot4,aes(x=plot4$Rank,y=plot4$percentage))+geom_point(size=3,color="blue")+
ylim(c(0, 1)) + geom_line(size=1.2,color="blue")+labs(x="Rank",y="percentage")
ggsave(p1,file="/home/qzzh/plot/bonemarrowandKidneyExp-barplot.pdf")
ggsave(p1,file="/home/qzzh/plot/bonemarrowandKidneyExp-barplot.ps")
ggsave(p2,file="/home/qzzh/plot/bonemarrowandKidneyExp-lineplot.pdf")
ggsave(p2,file="/home/qzzh/plot/bonemarrowandKidneyExp-lineplot.ps")
ggsave(p3,file="/home/qzzh/plot/bonemarrowandKidneyPhy-barplot.pdf")
ggsave(p3,file="/home/qzzh/plot/bonemarrowandKidneyPhy-barplot.ps")
ggsave(p4,file="/home/qzzh/plot/bonemarrowandKidneyPhy-lineplot.pdf")
ggsave(p4,file="/home/qzzh/plot/bonemarrowandKidneyPhy-lineplot.ps")
##################################################################################
##############################Figure6-split#######################################
##################################################################################
plot5 <- read.table("/home/qzzh/plot/neuron-and-nonneuronExp-barplot.txt",header = T,sep = "\t")
plot6 <- read.table("/home/qzzh/plot/neuron-and-nonneuronExp-lineplot.txt",header = T)
plot7 <- read.table("/home/qzzh/plot/neuron-and-nonneuronPhy-barplot.txt",header = T,sep = "\t")
plot8 <- read.table("/home/qzzh/plot/neuron-and-nonneuronPhy-lineplot.txt",header = T)
plot5 <- plot5[order(plot5$negativelog10),]
plot7 <- plot7[order(plot7$negativelog10),]
plot5$Desc <- with(plot5, factor(Desc, levels=Desc, ordered=TRUE))
plot7$Desc <- with(plot7, factor(Desc, levels=Desc, ordered=TRUE))
p5 <- ggplot(plot5,aes(x=plot5$Desc,y=plot5$negativelog10))+theme(panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(size=0.5),legend.title=element_blank())+geom_bar(fill="blue",stat="identity")+theme(axis.text.x=element_text(angle=0,colour="black",size=15,hjust=0.5))+theme(axis.text.y=element_text(angle=0,colour="black",size=10,hjust=1))+coord_flip()+labs(x="Desc",y="negativelog10(p-value)")+scale_y_discrete(limits=c(0,50,100,150,200))
p6 <- ggplot(plot6,aes(x=plot6$Rank,y=plot6$percentage))+geom_point(size=3,color="blue")+ylim(c(0, 1)) + geom_line(size=1.2,color="blue")+labs(x="Rank",y="percentage")
p7 <- ggplot(plot7,aes(x=plot7$Desc,y=plot7$negativelog10))+theme(panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(size=0.5),legend.title=element_blank())+geom_bar(fill="blue",stat="identity")+theme(axis.text.x=element_text(angle=0,colour="black",size=15,hjust=0.5))+theme(axis.text.y=element_text(angle=0,colour="black",size=10,hjust=1))+coord_flip()+labs(x="Desc",y="negativelog10(p-value)")+scale_y_discrete(limits=c(0,50,100,150,200))
p8 <- ggplot(plot8,aes(x=plot8$Rank,y=plot8$percentage))+geom_point(size=3,color="blue")+ylim(c(0,1)) + geom_line(size=1.2,color="blue")+labs(x="Rank",y="percentage")
ggsave(p5,file="/home/qzzh/plot/neuron-and-nonneuronExp-barplot.pdf")
ggsave(p5,file="/home/qzzh/plot/neuron-and-nonneuronExp-barplot.ps") 
ggsave(p6,file="/home/qzzh/plot/neuron-and-nonneuronExp-lineplot.pdf")
ggsave(p6,file="/home/qzzh/plot/neuron-and-nonneuronExp-lineplot.ps")
ggsave(p7,file="/home/qzzh/plot/neuron-and-nonneuronPhy-barplot.pdf")
ggsave(p7,file="/home/qzzh/plot/neuron-and-nonneuronPhy-barplot.ps") 
ggsave(p8,file="/home/qzzh/plot/neuron-and-nonneuronPhy-lineplot.pdf")
ggsave(p8,file="/home/qzzh/plot/neuron-and-nonneuronPhy-lineplot.ps")
