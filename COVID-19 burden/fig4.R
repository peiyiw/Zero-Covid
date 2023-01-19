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
library(RColorBrewer)

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
    axis.title.y.right =element_text(size = 12,color="black",vjust=0,angle=90),
    axis.title.x = element_text(size = 12, color="black",vjust = 0),
    axis.text.y.left = element_text(size = 12,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 12,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 12, color="black",vjust = 0.5, hjust=0.5, angle = 0),
    plot.margin=unit(rep(0.5,4),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
cl1 = rgb(99, 99, 99, maxColorValue = 255)

Sys.setlocale(category = "LC_CTYPE", locale = "chs")

filepath <- 'G:\\2023\\20230104_zero\\figure\\fig4\\plot_0116'
setwd(filepath)
colorset = c("#433983", "#20908c", "#34b778","#fde625")
capacity_hos <- 3080364
capacity_icu <- 63527

disname <- c('asthma','cancer','ckd','cld','copd','diabetes','hypertension','obesity')
dislevel <- c('hypertension','diabetes','ckd','obesity','copd','asthma','cld','cancer')
dislabel <- c('Hypertension','Diabetes','Chronic kidney disease','Obesity',
              'Chronic obstructive\npulmonary disease','Asthma','Chronic lung disease','Cancer')
hos_dataset <- data.frame()
icu_dataset <- data.frame()

interval <- 3
lag <- 21

for (i in c(1:length(disname))){
  # hospitalization
  hos_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_40_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_40_60_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_60_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  
  hos_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_40_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_40_60_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_60_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  
  hos_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_20_40_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_40_60_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  hos_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_dis_extra_60_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  
  
  hos_mean <- data.frame("20"=hos_mean_1[365],"20_40"=hos_mean_2[365],"40_60"=hos_mean_3[365],"60"=hos_mean_4[365])
  hos_up <- data.frame("20"=hos_up_1[365],"20_40"=hos_up_2[365],"40_60"=hos_up_3[365],"60"=hos_up_4[365])
  hos_down <- data.frame("20"=hos_down_1[365],"20_40"=hos_down_2[365],"40_60"=hos_down_3[365],"60"=hos_down_4[365])
  
  hos_mean <- melt(hos_mean,variable.name="age",value.name = "mean")
  hos_up <- melt(hos_up,variable.name="age",value.name = "up")
  hos_down <- melt(hos_down,variable.name="age",value.name = "down")
  
  hos_mean$dis <- disname[i]
  hos_up$dis <- disname[i]
  hos_down$dis <- disname[i]
  
  hos <- join(hos_mean,hos_up)
  hos <- join(hos,hos_down)

  hos_dataset <- rbind(hos_dataset,hos)
  
  # icu
  icu_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_40_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_40_60_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_60_mean_",i,".xlsx",sep=""),sheet=as.character(lag)))

  icu_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_40_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_40_60_up_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_60_up_",i,".xlsx",sep=""),sheet=as.character(lag)))

  icu_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_20_40_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_40_60_down_",i,".xlsx",sep=""),sheet=as.character(lag)))
  icu_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_dis_extra_60_down_",i,".xlsx",sep=""),sheet=as.character(lag)))

  
  icu_mean <- data.frame("20"=icu_mean_1[365],"20_40"=icu_mean_2[365],"40_60"=icu_mean_3[365],"60"=icu_mean_4[365])
  icu_up <- data.frame("20"=icu_up_1[365],"20_40"=icu_up_2[365],"40_60"=icu_up_3[365],"60"=icu_up_4[365])
  icu_down <- data.frame("20"=icu_down_1[365],"20_40"=icu_down_2[365],"40_60"=icu_down_3[365],"60"=icu_down_4[365])
  
  icu_mean <- melt(icu_mean,variable.name="age",value.name = "mean")
  icu_up <- melt(icu_up,variable.name="age",value.name = "up")
  icu_down <- melt(icu_down,variable.name="age",value.name = "down")
  
  icu_mean$dis <- disname[i]
  icu_up$dis <- disname[i]
  icu_down$dis <- disname[i]
  
  icu <- join(icu_mean,icu_up)
  icu <- join(icu,icu_down)
  
  icu_dataset <- rbind(icu_dataset,icu)
}

scale <- 100000

hos_dataset$dis<-factor(hos_dataset$dis,levels=rev(dislevel),labels = rev(dislabel))
icu_dataset$dis<-factor(icu_dataset$dis,levels=rev(dislevel),labels = rev(dislabel))

p_hos<-ggplot(data=hos_dataset)+
  geom_hline(yintercept = 0)+
  geom_bar(aes(x = dis, y = mean/scale, fill=age,group=age),stat = "identity",position=position_dodge())+
  geom_errorbar(aes(x = dis, ymax = up/scale, ymin = down/scale,group=age), 
                position = position_dodge(0.9),width=0.5,color=cl1,size=0.25)+
  coord_flip()+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x=NULL,y=paste("Extra Hospitalizations (100,000)",sep=""))+
  scale_y_continuous(breaks=breaks_pretty(4,min.n=4,max.n=4))+
  mytheme+theme(strip.text.x = element_blank())

p_icu<-ggplot(data=icu_dataset)+
  geom_hline(yintercept = 0)+
  geom_bar(aes(x = dis, y = mean/scale, fill=age,group=age),stat = "identity",position=position_dodge())+
  geom_errorbar(aes(x = dis, ymax = up/scale, ymin = down/scale,group=age),
                position = position_dodge(0.9),width=0.5,color=cl1,size=0.25)+
  coord_flip()+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x=NULL,y=paste("Extra Hospitalizations (100,000)",sep=""))+
  scale_y_continuous(breaks=breaks_pretty(4,min.n=4,max.n=4))+
  mytheme+theme(strip.text.x = element_blank())


p <- ggarrange(p_hos,p_icu,
               nrow=1,ncol=2,heights = c(4,1))
p

pdf(paste(filepath,"\\output\\figure4.pdf",sep=""),width=10,height=4)
print(p)
dev.off()

#calculation
# hospitalization
hos_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_mean.xlsx",sep=""),sheet=as.character(lag)))
hos_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_40_mean.xlsx",sep=""),sheet=as.character(lag)))
hos_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_40_60_mean.xlsx",sep=""),sheet=as.character(lag)))
hos_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_60_mean.xlsx",sep=""),sheet=as.character(lag)))

hos_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_up.xlsx",sep=""),sheet=as.character(lag)))
hos_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_40_up.xlsx",sep=""),sheet=as.character(lag)))
hos_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_40_60_up.xlsx",sep=""),sheet=as.character(lag)))
hos_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_60_up.xlsx",sep=""),sheet=as.character(lag)))

hos_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_down.xlsx",sep=""),sheet=as.character(lag)))
hos_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_20_40_down.xlsx",sep=""),sheet=as.character(lag)))
hos_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_40_60_down.xlsx",sep=""),sheet=as.character(lag)))
hos_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_hospitalization_60_down.xlsx",sep=""),sheet=as.character(lag)))


hos_mean <- data.frame("20"=hos_mean_1[365],"20_40"=hos_mean_2[365],"40_60"=hos_mean_3[365],"60"=hos_mean_4[365])
hos_up <- data.frame("20"=hos_up_1[365],"20_40"=hos_up_2[365],"40_60"=hos_up_3[365],"60"=hos_up_4[365])
hos_down <- data.frame("20"=hos_down_1[365],"20_40"=hos_down_2[365],"40_60"=hos_down_3[365],"60"=hos_down_4[365])

hos_mean <- melt(hos_mean,variable.name="age",value.name = "mean")
hos_up <- melt(hos_up,variable.name="age",value.name = "up")
hos_down <- melt(hos_down,variable.name="age",value.name = "down")

hos_mean$dis <- "Without underlying\nhealth condition"
hos_up$dis <- "Without underlying\nhealth condition"
hos_down$dis <- "Without underlying\nhealth condition"

hos <- join(hos_mean,hos_up)
hos <- join(hos,hos_down)


# icu
icu_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_mean.xlsx",sep=""),sheet=as.character(lag)))
icu_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_40_mean.xlsx",sep=""),sheet=as.character(lag)))
icu_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_40_60_mean.xlsx",sep=""),sheet=as.character(lag)))
icu_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_60_mean.xlsx",sep=""),sheet=as.character(lag)))

icu_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_up.xlsx",sep=""),sheet=as.character(lag)))
icu_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_40_up.xlsx",sep=""),sheet=as.character(lag)))
icu_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_40_60_up.xlsx",sep=""),sheet=as.character(lag)))
icu_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_60_up.xlsx",sep=""),sheet=as.character(lag)))

icu_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_down.xlsx",sep=""),sheet=as.character(lag)))
icu_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_20_40_down.xlsx",sep=""),sheet=as.character(lag)))
icu_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_40_60_down.xlsx",sep=""),sheet=as.character(lag)))
icu_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_cum_icu_60_down.xlsx",sep=""),sheet=as.character(lag)))


icu_mean <- data.frame("20"=icu_mean_1[365],"20_40"=icu_mean_2[365],"40_60"=icu_mean_3[365],"60"=icu_mean_4[365])
icu_up <- data.frame("20"=icu_up_1[365],"20_40"=icu_up_2[365],"40_60"=icu_up_3[365],"60"=icu_up_4[365])
icu_down <- data.frame("20"=icu_down_1[365],"20_40"=icu_down_2[365],"40_60"=icu_down_3[365],"60"=icu_down_4[365])

icu_mean <- melt(icu_mean,variable.name="age",value.name = "mean")
icu_up <- melt(icu_up,variable.name="age",value.name = "up")
icu_down <- melt(icu_down,variable.name="age",value.name = "down")

icu_mean$dis <- "Without underlying\nhealth condition"
icu_up$dis <- "Without underlying\nhealth condition"
icu_down$dis <- "Without underlying\nhealth condition"

icu <- join(icu_mean,icu_up)
icu <- join(icu,icu_down)

#total extra
round(sum(hos_dataset$mean)/1000000,2)
round(sum(hos_dataset$up)/1000000,2)
round(sum(hos_dataset$down)/1000000,2)

round(sum(icu_dataset$mean)/1000000,2)
round(sum(icu_dataset$up)/1000000,2)
round(sum(icu_dataset$down)/1000000,2)

#x-fold to 0-day-interval
round((sum(hos_dataset$mean)/1000000)/(sum(hos$mean)/1000000),2)
round((sum(hos_dataset$up)/1000000)/(sum(hos$down)/1000000),2)
round((sum(hos_dataset$down)/1000000)/(sum(hos$up)/1000000),2)

round((sum(icu_dataset$mean)/1000000)/(sum(icu$mean)/1000000),2)
round((sum(icu_dataset$up)/1000000)/(sum(icu$down)/1000000),2)
round((sum(icu_dataset$down)/1000000)/(sum(icu$up)/1000000),2)

# proportion
round(sum(hos_dataset[which(hos_dataset$age=="X40_60"),"mean"])/sum(hos_dataset$mean)*100,2)
round(sum(icu_dataset[which(icu_dataset$age=="X60"),"mean"])/sum(icu_dataset$mean)*100,2)

#total age extra
aggregate(mean~age,hos_dataset,FUN="sum")
aggregate(mean~age,icu_dataset,FUN="sum")
