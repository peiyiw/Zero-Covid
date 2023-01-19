rm(list = ls())
library(ggplot2) 
library(openxlsx)
library(lubridate)
library(reshape2)
library(scales)
library(Rmisc)
library(ggpubr)
library(dplyr)
library(ggforce)


Sys.setlocale("LC_TIME", "English") #print in english

mytheme <- theme_minimal()+
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    plot.title=element_text(size=12, vjust = 0.5, hjust = 0.5),
    legend.background = element_rect(fill=NA, size=0,color=NA),
    legend.text=element_text(size=12),
    legend.position="none",
    legend.title=element_text(face="bold",size=12),
    axis.line.y=element_line(linetype=1,color='black',size=0),
    axis.line.x=element_line(linetype=1,color='black',size=0),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=2,color='grey'),
    axis.title.y.left = element_text(size = 12,color="black",vjust=2),
    axis.title.y.right =element_text(size = 12,color="black",vjust=0.5,angle=90),
    axis.title.x = element_text(size = 12, color="black",vjust = 0),
    axis.text.y.left = element_text(size = 12,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 12,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 12, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"), axis.ticks.y=element_line(color="black",size=.5),
    axis.ticks.x=element_line(color="black",size=.5),
    strip.text.x = element_text(size=12, color="black"),
    strip.text.y = element_text(size=12, color="black"),
    plot.margin=unit(rep(0.5,4),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
c = rgb(179, 179, 179, maxColorValue = 255)
c = rgb(229, 216, 189, maxColorValue = 255)
cl1 = rgb(217, 217, 217, maxColorValue = 255)
filepath <- 'G:\\2023\\20230104_zero\\figure\\fig2\\plot_0116'
setwd(filepath)

interval <- c(1,2,3,4)
lag <- seq(7,21,7)
data <- data.frame()
for (i in interval){
  for (j in lag){
    ## duration
    dur_mean <- read.xlsx(paste("..\\data\\interval_",i,"_duration_mean.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    dur_mean <- melt(dur_mean,variable.name="city",value.name="mean")
    dur_mean$type <- "duration (days)"
    dur_mean$lag <- as.character(j)
    dur_mean$interval <- as.character(i)
    
    dur_up <- read.xlsx(paste("..\\data\\interval_",i,"_duration_up.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    dur_up <- melt(dur_up,variable.name="city",value.name="up")
    dur_up$lag <- as.character(j)
    dur_up$interval <- as.character(i)
    
    dur_down <- read.xlsx(paste("..\\data\\interval_",i,"_duration_down.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    dur_down <- melt(dur_down,variable.name="city",value.name="down")
    dur_down$lag <- as.character(j)
    dur_down$interval <- as.character(i)
    
    dur <- merge(dur_mean,dur_up)
    dur <- merge(dur,dur_down)
    ## round
    round_mean <- read.xlsx(paste("..\\data\\interval_",i,"_round_mean.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    round_mean <- melt(round_mean,variable.name="city",value.name="mean")
    round_mean$type <- "round"
    round_mean$lag <- as.character(j)
    round_mean$interval <- as.character(i)
    
    round_up <- read.xlsx(paste("..\\data\\interval_",i,"_round_up.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    round_up <- melt(round_up,variable.name="city",value.name="up")
    round_up$lag <- as.character(j)
    round_up$interval <- as.character(i)
    
    round_down <- read.xlsx(paste("..\\data\\interval_",i,"_round_down.xlsx",sep=""), sheet=as.character(j),rowNames = T)
    round_down <- melt(round_down,variable.name="city",value.name="down")
    round_down$lag <- as.character(j)
    round_down$interval <- as.character(i)
    
    round <- merge(round_mean,round_up)
    round <- merge(round,round_down)
    
    ##data and confidence interval
    data <- rbind(data,dur,round)
  }
}
data <- data[-which(data$mean==0),]

# data[which(data$type=="dur_meanation (days)"&data$value>180),"value"] <- 180
# data[which(data$type=="round_means"&data$value>60),"value"] <- 60

data$lag <- factor(data$lag,unique(data$lag))
data$interval <- factor(data$interval, levels=interval, labels = c(paste(interval, "-day-interval",sep="")))

counts <- data %>% count(interval,lag,type)
counts <- counts[which(counts$type=="duration (days)"),]
data_hline <- data.frame(type = unique(data$type),
                         hline = c(76,NA))
tmp <- aggregate(mean~lag+interval+type,data=data,FUN=min)
counts <- cbind(counts, data.frame("value"=tmp[which(tmp$type=="duration (days)"),"mean"]))

geom.text.size = 10
theme.size = (5/14) * geom.text.size

# data <- data[-which(data$value>100&data$type=="dur_meanation (days)"),]
# data <- data[-which(data$value>40&data$type=="round_means"),]

# p<-ggplot(data=data,mapping=aes(x=lag, y=mean))+
#   facet_grid(type ~ interval,scales = "free_y",space='free_x')+
#   geom_violin(trim = T,color="black",fill="#f2f2f2")+
#   geom_jitter(mapping=aes(color=lag),width = 0.05,alpha=1,size=0.5,shape=21,color="black",fill="grey",stroke=0.05)+
#   geom_errorbar(aes(x=lag, ymin = down, ymax = up),width = 0)+
#   geom_text(data=counts,size=theme.size,aes(x=lag,y=value-25,label=paste("n=",n,sep=""))) +
#   # geom_hline(data=data_hline,aes(yintercept = hline),linetype=2,color="#b2182b",size=0.8) + 
#   scale_x_discrete(breaks=c(lag), labels=c("1","2","3"),name="Response lag (weeks)")+
#   scale_y_continuous(breaks = scales::pretty_breaks(3),expand=expansion(mult = c(0.15, 0.15)))+
#   # scale_color_manual(values=c(c1,c2,c3,c4))+
#   # scale_color_manual(values=rep(c,4))+
#   mytheme

p<-ggplot(data=data,mapping=aes(x=lag, y=mean))+
  facet_grid(type ~ interval,scales = "free_y",space='free_x')+
  geom_pointrange(mapping=aes(ymin = down, ymax = up,color=lag),
                  position=position_jitter(width = 0.15),
                  linetype="solid",
                  alpha=1,size=0.2,shape=21,
                  color=cl1,fill="darkgrey",stroke=0.03)+
  geom_violin(trim = T,color="black",fill=NA)+
  # geom_errorbar(aes(x=lag, ymin = down, ymax = up),width = 0)+
  geom_text(data=counts,size=theme.size,aes(x=lag,y=value-25,label=paste("n=",n,sep=""))) +
  geom_hline(data=data_hline,aes(yintercept = hline),linetype=2,color="#b2182b",size=0.8) +
  scale_x_discrete(breaks=c(lag), labels=c("1","2","3"),name="Response lag (weeks)")+
  scale_y_continuous(breaks = scales::pretty_breaks(4),expand=expansion(mult = c(0.15, 0.15)))+
  # scale_color_manual(values=c(c1,c2,c3,c4))+
  # scale_color_manual(values=rep(c,4))+
  mytheme

p

pdf(paste(filepath,"\\output\\fig3_a.pdf",sep=""),width=5.5,height=3.42)
print(p)
dev.off()

