rm(list = ls())
############
# Rudolf Cardinal, March 2011
# Simple extensions to ggplot2 (v0.8.7); see http://pobox.com/~rudolf/statistics/R
# Modified 5 Jan 2013 for ggplot2 0.9.3 (NB: use sessionInfo() to find current package versions)
# - fetch ggplot2 source with: git clone https://github.com/hadley/ggplot2.git
# Changes, because ggplot2 has changed its internal calling mechanisms:
# - opts() deprecated in favour of theme()
# - "Element panel.border must be an element_rect object" (error from validate_element() in theme-elements.r)
#   ... so change all class = "theme" to class = c("element_rect", "element")
# - "cannot coerce type 'closure' to vector of type 'list'"
#   ... a closure is a function (see ?typeof)
#   ... change class to be of class c("MYCLASS", "element_rect", "element")
# - then element_grob.MYCLASS not called by element_render()/element_grob()/UseMethod()... environment/namespace problem
#   tried setMethod("element_grob", "theme_border", function(STUFF) { STUFF} , where = as.environment("package:ggplot2")
#   but the environment is locked
#   ggplot2's theme-elements.r defines e.g. element_rect (exported) and element_grob.element_rect (not exported, does the work)
#   However, we can't override an internal function:
#       ... e.g. rewrite "validate_element" to crash
#           set environment(validate_element) <- as.environment("package:ggplot2") -- doesn't break the plotting.
# - Upshot: now impossible to hack through this way (locked environment).
# - http://obeautifulcode.com/R/How-R-Searches-And-Finds-Stuff/
# - http://stackoverflow.com/questions/8204008/redirect-intercept-function-calls-within-a-package-function
# - These don't fix it:
#   library(proto)
#   theme <- with(proto(environment(ggplot2::theme), theme = ggplot2::theme, element_grob.theme_border = my.element_grob.theme_border), theme) --- doesn't work
#   ggplot <- with(proto(environment(ggplot2::ggplot), ggplot = ggplot2::ggplot, element_grob.theme_border = my.element_grob.theme_border), ggplot) --- breaks!
# - Fix by Baptiste Auguie 8/1/2013: inherit from element_blank instead; then it works fine.

#-------------------------------------------------------------------------------
# Requirements
#-------------------------------------------------------------------------------

library(grid) # for gpar

#-------------------------------------------------------------------------------
# Code duplicated from ggplot2 source (not exposed to wider namespace) for convenience
#-------------------------------------------------------------------------------

.pt <- 1 / 0.352777778
len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}

#-------------------------------------------------------------------------------
# Generic panel border (can set any combination of left/right/top/bottom)
#-------------------------------------------------------------------------------

theme_border <- function(
  type = c("left", "right", "bottom", "top", "none"),
  colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + opts( panel.border=theme_border(type=c("bottom","left")) ) + ...
  type <- match.arg(type, several.ok=TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}
element_grob.theme_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (length(type)==0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: "L" (left + bottom) border
#-------------------------------------------------------------------------------

theme_L_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_L_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_L_border", "element_blank", "element")
  )
}
element_grob.theme_L_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x+width, x, x), y = c(y,y,y+height), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: bottom border only
#-------------------------------------------------------------------------------

theme_bottom_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_bottom_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_bottom_border", "element_blank", "element")
  )
}
element_grob.theme_bottom_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x, x+width), y = c(y,y), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: left border only
#-------------------------------------------------------------------------------

theme_left_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_left_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_left_border", "element_blank", "element")
  )
}
element_grob.theme_left_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x, x), y = c(y, y+height), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}



#-------------------------------------------------------------------------------
# Border selection by number
#-------------------------------------------------------------------------------
theme_border_numerictype <- function(type, colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_border(type=9) ) + ...
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border_numerictype", "element_blank", "element")
  )
}
element_grob.theme_border_numerictype <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  # numerical types from: library(gridExtra); example(borderGrob)
  # 1=none, 2=bottom, 3=right, 4=top, 5=left, 6=B+R, 7=T+R, 8=T+L, 9=B+L, 10=T+B, 11=L+R, 12=T+B+R, 13=T+L+R, 14=T+B+L, 15=B+L+R, 16=T+B+L+R
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if (type==2 || type==6 || type==9 || type==10 || type==12 || type==14 || type==15 || type==16) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if (type==4 || type==7 || type==8 || type==10 || type==12 || type==13 || type==14 || type==16) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if (type==5 || type==8 || type==9 || type==11 || type==13 || type==14 || type==15 || type==16) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if (type==3 || type==6 || type==7 || type==11 || type==12 || type==13 || type==15 || type==16) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (type==1) { # blank; can't pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# Examples
#-------------------------------------------------------------------------------

rnc_ggplot2_border_themes_example_script = '
    library(ggplot2)
    df = data.frame( x=c(1,2,3), y=c(4,5,6) )
    source("http://egret.psychol.cam.ac.uk/statistics/R/extensions/rnc_ggplot2_border_themes_2013_01.r")
    ggplot(data=df, aes(x=x, y=y)) + geom_point() + theme_bw() + theme( panel.border = theme_border( c("bottom","left") ) )
    ggplot(data=df, aes(x=x, y=y)) + geom_point() + theme_bw() + theme( panel.border = theme_left_border() )
    ggplot(data=df, aes(x=x, y=y)) + geom_point() + theme_bw() + theme( panel.border = theme_bottom_border() )
    ggplot(data=df, aes(x=x, y=y)) + geom_point() + theme_bw() + theme( panel.border = theme_L_border() )
    ggplot(data=df, aes(x=x, y=y)) + geom_point() + theme_bw() + theme( panel.border = theme_border_numerictype(12) ) # use 1:16 as possible values
'

#######

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
    axis.line.y=element_blank(),
    axis.line.x=element_blank(),
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
    plot.margin=unit(c(0.25,0.25,0.25,0.25),'lines')
  )


filepath <- 'G:\\2022\\20220907_testing_adjust\\figure\\fig3\\plot_1021'
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
  for (interval in c(0,3)){
    svname <- outcome[isv]
    #####
    ## read data
    sv1 <- read.xlsx(paste("..\\data\\interval_",interval,"_",filename[isv],"_20.xlsx",sep=""), sheet=as.character(lag))
    sv2 <- read.xlsx(paste("..\\data\\interval_",interval,"_",filename[isv],"_20_40.xlsx",sep=""), sheet=as.character(lag))
    sv3 <- read.xlsx(paste("..\\data\\interval_",interval,"_",filename[isv],"_40_60.xlsx",sep=""), sheet=as.character(lag))
    sv4 <- read.xlsx(paste("..\\data\\interval_",interval,"_",filename[isv],"_60.xlsx",sep=""), sheet=as.character(lag))
    city_attri <- read.xlsx("..\\data\\city_attribute.xlsx", sheet="city")
    nation_age_structure <- read.xlsx("..\\data\\data_nation_age_structure.xlsx", sheet="data")
    pro_age_structure <- read.xlsx("..\\data\\data_pro_age_structure.xlsx", sheet="data")
    
    sv <- sv1 + sv2 + sv3 + sv4
    sv <- rowSums(sv)
    peakday <- which(sv==max(sv))
    # plot(sv)
    
    sv1$day <- c(1:dim(sv1)[1])
    sv2$day <- c(1:dim(sv2)[1])
    sv3$day <- c(1:dim(sv3)[1])
    sv4$day <- c(1:dim(sv4)[1])
    
    sv1 <- melt(sv1,id.vars = "day",variable.name = "name",value.name = "baseline20")
    sv2 <- melt(sv2,id.vars = "day",variable.name = "name",value.name = "baseline20_40")
    sv3 <- melt(sv3,id.vars = "day",variable.name = "name",value.name = "baseline40_60")
    sv4 <- melt(sv4,id.vars = "day",variable.name = "name",value.name = "baseline60")
    
    sv_age <- merge(sv1,sv2)
    sv_age <- merge(sv_age,sv3)
    sv_age <- merge(sv_age,sv4)
    
    sv_age_city <- merge(city_attri,sv_age)
    baseline_sv_dataset <- aggregate(cbind(baseline20,baseline20_40,baseline40_60,baseline60)~(day+provincech+provinceen),data=sv_age_city,sum)
    
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
        dis_sv <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
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
        dis_sv <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
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
        dis_sv <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
        dis_sv_extra <- data.frame("day"=tmp[,"day"],"provinceen"=tmp[,"provinceen"])
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
    
    nation_baseline_sv <- baseline_sv_dataset
    nation_dis_sv_extra <- dis_sv_extra_dataset
    
    nation_baseline_sv <- aggregate(cbind(baseline20,baseline20_40,baseline40_60,baseline60)~(day),nation_baseline_sv,sum)
    nation_dis_sv_extra <- aggregate(cbind(extra20,extra20_40,extra40_60,extra60)~(day),nation_dis_sv_extra,sum)
    
    tmp <- merge(nation_dis_sv_extra,nation_baseline_sv)
    nation_dis_sv <- data.frame("day"=tmp$day,
                                      "dis20"=tmp$extra20+tmp$baseline20,
                                      "dis20_40"=tmp$extra20_40+tmp$baseline20_40,
                                      "dis40_60"=tmp$extra40_60+tmp$baseline40_60,
                                      "dis60"=tmp$extra60+tmp$baseline60)
    
    nation_baseline_sv <- melt(nation_baseline_sv,id.vars=c("day"),variable="age",value.name = paste(svname,sep=""))
    nation_dis_sv <- melt(nation_dis_sv,id.vars=c("day"),variable="age",value.name = paste(svname,sep=""))
    
    assign(paste("peakday_",svname,"_interval_",interval,sep=""),peakday)
    assign(paste("nation_baseline_",svname,"_interval_",interval,sep=""),nation_baseline_sv)
    assign(paste("nation_dis_",svname,"_interval_",interval,sep=""),nation_dis_sv)
  }
}

nation_dis_hos_interval_0 <- merge(nation_dis_hos_interval_0,nation_dis_icu_interval_0)
nation_dis_hos_interval_0$hos_all <- nation_dis_hos_interval_0$hos+
  nation_dis_hos_interval_0$icu
nation_dis_hos_interval_3 <- merge(nation_dis_hos_interval_3,nation_dis_icu_interval_3)
nation_dis_hos_interval_3$hos_all <- nation_dis_hos_interval_3$hos+
  nation_dis_hos_interval_3$icu


scale <- 1000

p1 <- ggplot()+
  geom_bar(data=nation_dis_hos_interval_3,aes(x = day, y = hos_all/scale, fill=age,group=age), stat = "identity")+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x="Days since arrival of COVID-19 in China",y=paste("Hospital beds required (1,000)",sep=""))+
  scale_y_continuous(limits=c(0,220),breaks=c(0,100,200),labels = c("0","100","     200"),
                     sec.axis = sec_axis(~.*scale/capacity_hos*100,name = "Hospital occupancy (%)",
                                         breaks=c(0,3,6),
                                         labels = c("0","3","6        ")))+
  mytheme+theme(panel.border = theme_border(type = c("right","bottom","left")))
p1

p2 <- ggplot()+
  geom_hline(yintercept=capacity_hos/scale,
             color="#b2182b", linetype="dashed",size=0.8)+
  geom_hline(yintercept=max(sum(nation_dis_hos_interval_0[which(nation_dis_hos_interval_0$day==peakday_hos_interval_0),"hos_all"]),
             sum(nation_dis_hos_interval_0[which(nation_dis_hos_interval_0$day==peakday_icu_interval_0),"hos_all"]))/scale,
             color="grey", linetype="dashed",size=0.8)+
  scale_y_continuous(limits=c(2800,12000),breaks=c(3000,11000),
                     sec.axis = sec_axis(~.*scale/capacity_hos*100,
                                         name = "Hospital occupancy (%)",breaks=c(0,100,350),labels = c("0","100","350    ")))+
  mytheme+theme(panel.border = theme_border(type = c("top","left","right")))
p2


p_hos <- ggarrange(p2,p1,nrow=2,heights = c(1,3))
p_hos


p3 <- ggplot(data=nation_dis_icu_interval_3)+
  geom_bar(aes(x = day, y = icu/scale, fill=age,group=age), stat = "identity")+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x="Days since arrival of COVID-19 in China",y=paste("ICU beds required (1,000)",sep=""))+
  scale_y_continuous(limits=c(0,40),breaks=c(0,18,36), labels = c("0","20","    40"),
                     sec.axis = sec_axis(~.*scale/capacity_icu*100,name = "ICU occupancy (%)",
                                         breaks=c(0,25,50),labels = c("0","25","50    ")))+
  mytheme+theme(panel.border = theme_border(type = c("right","bottom","left")))
p3



p4 <- ggplot()+
  geom_hline(yintercept=capacity_icu/scale,
             color="#b2182b", linetype="dashed",size=0.8)+
  geom_hline(yintercept=sum(nation_dis_icu_interval_0[which(nation_dis_icu_interval_0$day==peakday_icu_interval_0),"icu"])/scale,
             color="grey", linetype="dashed",size=0.8)+
  scale_y_continuous(limits=c(20,1700), breaks=c(100,1600), labels = c("100","1600"),
                     sec.axis = sec_axis(~.*scale/capacity_icu*100,name = "ICU occupancy (%)",breaks=c(100,2500)))+
  mytheme+theme(panel.border = theme_border(type = c("top","left","right")))
p4


p_icu <- ggarrange(p4,p3,nrow=2,heights = c(1,3))
p_icu

p <- ggarrange(p_hos,p_icu,ncol=2)
p

pdf(paste(filepath,"\\output\\fig3_a.pdf",sep=""),width=10,height=3.5)
print(p)
dev.off()

