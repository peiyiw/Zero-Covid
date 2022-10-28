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


Sys.setlocale("LC_TIME", "English") #print in english

#####
getPalette = colorRampPalette(brewer.pal(6, "YlGnBu"))
mytheme <- theme_minimal()+
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    plot.title=element_text(size=12, vjust = 0.5, hjust = 0),
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
    strip.text.y = element_text(size=12, color="black", hjust = 0),
    plot.margin=unit(rep(0.5,4),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

filepath <- 'G:\\2022\\20220907_testing_adjust\\figure\\fig4\\plot_1021'
setwd(filepath)


## data
ngroup <- 4
lag=21
hos_or <- data.frame("obesity"=5.82,"asthma"=3.83,"cancer"=2.88,
                     "CKD"=3.4,"COPD"=1.07,"diabetes"=2.27,
                     "CLD"=2.77,"hypertension"=2.01)
icu_or <- data.frame("obesity"=1.43,"asthma"=2.63,"cancer"=0.91,
                     "CKD"=5.32,"COPD"=6.66,"diabetes"=3.07,
                     "CLD"=1.1,"hypertension"=2.95)
or <- rbind(hos_or,icu_or)
rownames(or) <- c("hos","icu")
filename <- c("hospitalization","icu")
outcome <- c("hos","icu")
capacity_hos <- 3080364
capacity_icu <- 63527


for (isv in c(1:length(outcome))){
    svname <- outcome[isv]
    #####
    ## read data
    sv1 <- read.xlsx(paste("..\\data\\interval_3_cum_",filename[isv],"_20.xlsx",sep=""), sheet=as.character(lag))
    sv2 <- read.xlsx(paste("..\\data\\interval_3_cum_",filename[isv],"_20_40.xlsx",sep=""), sheet=as.character(lag))
    sv3 <- read.xlsx(paste("..\\data\\interval_3_cum_",filename[isv],"_40_60.xlsx",sep=""), sheet=as.character(lag))
    sv4 <- read.xlsx(paste("..\\data\\interval_3_cum_",filename[isv],"_60.xlsx",sep=""), sheet=as.character(lag))
    city_attri <- read.xlsx("..\\data\\city_attribute.xlsx", sheet="city")
    nation_age_structure <- read.xlsx("..\\data\\data_nation_age_structure.xlsx", sheet="data")
    pro_age_structure <- read.xlsx("..\\data\\data_pro_age_structure.xlsx", sheet="data")
    
    
    sv_age <- data.frame(t(sv1),t(sv2),t(sv3),t(sv4))
    colnames(sv_age) <- c("baseline20","baseline20_40","baseline40_60","baseline60")
    sv_age$name <- rownames(sv_age)
    
    sv_age_city <- merge(city_attri,sv_age)
    baseline_sv_dataset <- aggregate(cbind(baseline20,baseline20_40,baseline40_60,baseline60)~(provincech+provinceen),data=sv_age_city,sum)
    
    disease <- c("obesity", "diabetes","CKD","hypertension")
    agegroup <- c("20","20_40","40_60","60")
    
    dis_sv_dataset <- data.frame()
    dis_sv_extra_dataset <- data.frame()
    #####
    ## clinical severity
    for (dis in disease){
      dis_sv_or <- or[svname,dis]
      #read prevalence
      if (dis %in% c("cancer","COPD"))  { ##only pro
        dis_pr_pro <- read.xlsx(paste("..\\data\\",dis,"_data_pro.xlsx",sep=""), sheet="data")
        
        #case with disease
        tmp <- merge(dis_pr_pro,baseline_sv_dataset)
        dis_sv <- data.frame("provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("provinceen"=tmp[,"provinceen"])
        for (group in agegroup){
          dis_sv[,paste("dis",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            tmp[,"pr"]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-tmp[,"pr"])
          dis_sv_extra[,paste("extra",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            tmp[,"pr"]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-tmp[,"pr"])-
            tmp[,paste("baseline",group,sep="")]
        }
      }  else if (dis == "CKD"){ ##only age
        dis_pr_age <- read.xlsx(paste("..\\data\\",dis,"_data_age.xlsx",sep=""), sheet="data")
        
        #case with disease
        tmp <- baseline_sv_dataset
        dis_sv <- data.frame("provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("provinceen"=tmp[,"provinceen"])
        for (group in agegroup){
          dis_sv[,paste("dis",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            dis_pr_age[,group]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-dis_pr_age[,group])
          dis_sv_extra[,paste("extra",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            dis_pr_age[,group]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-dis_pr_age[,group])-
            tmp[,paste("baseline",group,sep="")]
        }
      }  else {
        dis_pr_age <- read.xlsx(paste("..\\data\\",dis,"_data_age.xlsx",sep=""), sheet="data")
        dis_pr_pro <- read.xlsx(paste("..\\data\\",dis,"_data_pro.xlsx",sep=""), sheet="data")
        age_num <- dis_pr_age*nation_age_structure
        age_num <- age_num/sum(age_num)
        
        tmp <- merge(dis_pr_pro,pro_age_structure)
        dis_pr_age_pro <- data.frame("provincech"=tmp[,"provincech"])
        for (group in agegroup){
          dis_pr_age_pro[,paste("pr",group,sep="")] <- pmin(tmp$pr*age_num[,group]/tmp[,group],1)
        }
        
        #case with disease
        tmp <- merge(dis_pr_age_pro,baseline_sv_dataset)
        dis_sv <- data.frame("provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("provinceen"=tmp[,"provinceen"])
        for (group in agegroup){
          dis_sv[,paste("dis",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            tmp[,paste("pr",group,sep="")]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-tmp[,paste("pr",group,sep="")])
          dis_sv_extra[,paste("extra",group,sep="")] <- tmp[,paste("baseline",group,sep="")]*
            tmp[,paste("pr",group,sep="")]*dis_sv_or+tmp[,paste("baseline",group,sep="")]*(1-tmp[,paste("pr",group,sep="")])-
            tmp[,paste("baseline",group,sep="")]
        }
      }
      
      dis_sv$type <- dis
      dis_sv_extra$type <- dis
      
      dis_sv_dataset <- rbind(dis_sv_dataset,dis_sv)
      dis_sv_extra_dataset <- rbind(dis_sv_extra_dataset,dis_sv_extra)
    }
    
    nation_dis_sv_extra <- dis_sv_extra_dataset
    nation_dis_sv_extra <- aggregate(cbind(extra20,extra20_40,extra40_60,extra60)~(type),nation_dis_sv_extra,sum)
    nation_dis_sv_extra <- melt(nation_dis_sv_extra,id.vars=c("type"),variable="age",value.name = paste(svname,sep=""))
    nation_dis_sv_extra$type <- factor(nation_dis_sv_extra$type,levels=c("obesity", "CKD","diabetes","hypertension"))
    
    assign(paste("nation_dis_",svname,"_extra",sep=""),nation_dis_sv_extra)
}

nation_dis_hos_extra <- merge(nation_dis_hos_extra,nation_dis_icu_extra)
nation_dis_hos_extra$hos_all <- nation_dis_hos_extra$hos+
  nation_dis_hos_extra$icu

scale <- 1000
###  plot
# extra
p_hos <- ggplot(data=nation_dis_hos_extra,aes(x = type, y = hos_all/scale, fill=age,group=age))+
  geom_hline(yintercept = 0)+
  geom_bar(stat = "identity",position=position_dodge())+
  coord_flip()+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x=NULL,y=paste("Extra Hospitalizations (1,000)",sep=""))+
  scale_y_continuous(breaks=breaks_pretty(4,min.n=4,max.n=4))+
  mytheme+theme(strip.text.x = element_blank())

p_icu <- ggplot(data=nation_dis_icu_extra,aes(x = type, y = icu/scale, fill=age,group=age))+
  geom_hline(yintercept = 0)+
  geom_bar(stat = "identity",position=position_dodge())+
  coord_flip()+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x=NULL,y=paste("Extra ICU admissions (1,000)",sep=""))+
  scale_y_continuous(breaks=breaks_pretty(4,min.n=4,max.n=4))+
  mytheme+theme(strip.text.x = element_blank())

p <- ggarrange(p_hos,p_icu,
          nrow=1,ncol=2,heights = c(4,1))
p

pdf(paste(filepath,"\\output\\figure4a.pdf",sep=""),width=12,height=4)
print(p)
dev.off()

