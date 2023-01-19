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
library(stringr)

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
Sys.setlocale("LC_TIME", "English") #print in english

filepath <- 'G:\\2023\\20230104_zero\\figure\\fig3\\plot_0119'
setwd(filepath)

colorset = c("#433983", "#20908c", "#34b778","#fde625")
capacity_hos <- 3080364
capacity_icu <- 13.81*10000

lag=21
cl1 = rgb(99, 99, 99, maxColorValue = 255)
cl2 = rgb(217, 217, 217, maxColorValue = 255)
scale <- 100000
## baseline data
for (interval in c(0,3,4)){
# for (interval in c(3)){
  if (interval==0){lag = 7}
  if (interval==3){lag = 21}
  # hospitalization
  hos_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_mean.xlsx",sep=""), sheet=as.character(lag)))
  hos_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_40_mean.xlsx",sep=""), sheet=as.character(lag)))
  hos_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_40_60_mean.xlsx",sep=""), sheet=as.character(lag)))
  hos_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_60_mean.xlsx",sep=""), sheet=as.character(lag)))

  hos_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_up.xlsx",sep=""), sheet=as.character(lag)))
  hos_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_40_up.xlsx",sep=""), sheet=as.character(lag)))
  hos_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_40_60_up.xlsx",sep=""), sheet=as.character(lag)))
  hos_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_60_up.xlsx",sep=""), sheet=as.character(lag)))
  
  hos_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_down.xlsx",sep=""), sheet=as.character(lag)))
  hos_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_20_40_down.xlsx",sep=""), sheet=as.character(lag)))
  hos_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_40_60_down.xlsx",sep=""), sheet=as.character(lag)))
  hos_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_hospitalization_dis_60_down.xlsx",sep=""), sheet=as.character(lag)))
  
  
  hos_mean <- data.frame("date"=c(1:length(hos_mean_1)),"20"=hos_mean_1,"20_40"=hos_mean_2,
                     "40_60"=hos_mean_3,"60"=hos_mean_4)
  hos_up <- data.frame("date"=c(1:length(hos_up_1)),"20"=hos_up_1,"20_40"=hos_up_2,
                         "40_60"=hos_up_3,"60"=hos_up_4)
  hos_down <- data.frame("date"=c(1:length(hos_down_1)),"20"=hos_down_1,"20_40"=hos_down_2,
                         "40_60"=hos_down_3,"60"=hos_down_4)
  
  
  hos_mean <- melt(hos_mean,id.vars = "date",variable.name="age",value.name = "mean")
  hos_up <- melt(hos_up,id.vars = "date",variable.name="age",value.name = "up")
  hos_down <- melt(hos_down,id.vars = "date",variable.name="age",value.name = "down")
  
  hos <- join(hos_mean,hos_up)
  hos <- join(hos,hos_down)
  
  if (interval==3){
    hos <- hos[-which(hos$up<1&hos$date>20),]
  }
  hos_ci <- aggregate(cbind(mean,up,down)~date,hos,sum)
  
  peakday_mean <- hos[which(hos$mean==max(hos$mean)),"date"]
  peakday_up <- hos[which(hos$up==max(hos$up)),"date"]
  peakday_down <- hos[which(hos$down==max(hos$down)),"date"]
  assign(paste("hos_peakday_mean_interval_",interval,sep=""),peakday_mean)
  assign(paste("hos_peakday_up_interval_",interval,sep=""),peakday_up)
  assign(paste("hos_peakday_down_interval_",interval,sep=""),peakday_down)
  assign(paste("hos_interval_",interval,sep=""),hos)
  assign(paste("hos_ci_interval_",interval,sep=""),hos_ci)
  
  # icu
  icu_mean_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_mean.xlsx",sep=""), sheet=as.character(lag)))
  icu_mean_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_40_mean.xlsx",sep=""), sheet=as.character(lag)))
  icu_mean_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_40_60_mean.xlsx",sep=""), sheet=as.character(lag)))
  icu_mean_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_60_mean.xlsx",sep=""), sheet=as.character(lag)))

  icu_up_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_up.xlsx",sep=""), sheet=as.character(lag)))
  icu_up_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_40_up.xlsx",sep=""), sheet=as.character(lag)))
  icu_up_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_40_60_up.xlsx",sep=""), sheet=as.character(lag)))
  icu_up_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_60_up.xlsx",sep=""), sheet=as.character(lag)))
  
  icu_down_1 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_down.xlsx",sep=""), sheet=as.character(lag)))
  icu_down_2 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_20_40_down.xlsx",sep=""), sheet=as.character(lag)))
  icu_down_3 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_40_60_down.xlsx",sep=""), sheet=as.character(lag)))
  icu_down_4 <- rowSums(read.xlsx(paste("..\\data\\interval_",interval,"_icu_dis_60_down.xlsx",sep=""), sheet=as.character(lag)))
  
  
  icu_mean <- data.frame("date"=c(1:length(icu_mean_1)),"20"=icu_mean_1,"20_40"=icu_mean_2,
                         "40_60"=icu_mean_3,"60"=icu_mean_4)
  icu_up <- data.frame("date"=c(1:length(icu_up_1)),"20"=icu_up_1,"20_40"=icu_up_2,
                       "40_60"=icu_up_3,"60"=icu_up_4)
  icu_down <- data.frame("date"=c(1:length(icu_down_1)),"20"=icu_down_1,"20_40"=icu_down_2,
                         "40_60"=icu_down_3,"60"=icu_down_4)
  
  
  icu_mean <- melt(icu_mean,id.vars = "date",variable.name="age",value.name = "mean")
  icu_up <- melt(icu_up,id.vars = "date",variable.name="age",value.name = "up")
  icu_down <- melt(icu_down,id.vars = "date",variable.name="age",value.name = "down")
  
  icu <- join(icu_mean,icu_up)
  icu <- join(icu,icu_down)
  
  if (interval==3){
    icu <- icu[-which(icu$up<1&icu$date>20),]
  }
  icu_ci <- aggregate(cbind(mean,up,down)~date,icu,sum)
  
  peakday_mean <- icu[which(icu$mean==max(icu$mean)),"date"]
  peakday_up <- icu[which(icu$up==max(icu$up)),"date"]
  peakday_down <- icu[which(icu$down==max(icu$down)),"date"]
  assign(paste("icu_peakday_mean_interval_",interval,sep=""),peakday_mean)
  assign(paste("icu_peakday_up_interval_",interval,sep=""),peakday_up)
  assign(paste("icu_peakday_down_interval_",interval,sep=""),peakday_down)
  assign(paste("icu_interval_",interval,sep=""),icu)
  assign(paste("icu_ci_interval_",interval,sep=""),icu_ci)
}


# hospitalization
xmin <- 0
xmax <- max(hos_interval_3$date)
hos_interval_4_mean_peak <- sum(hos_interval_4[which(hos_interval_4$date==hos_peakday_mean_interval_4),"mean"])
hos_interval_4_up_peak <- sum(hos_interval_4[which(hos_interval_4$date==hos_peakday_mean_interval_4),"up"])
hos_interval_4_down_peak <- sum(hos_interval_4[which(hos_interval_4$date==hos_peakday_mean_interval_4),"down"])
hos_interval_0_mean_peak <- sum(hos_interval_0[which(hos_interval_0$date==hos_peakday_mean_interval_0),"mean"])
hos_interval_0_up_peak <- sum(hos_interval_0[which(hos_interval_0$date==hos_peakday_mean_interval_0),"up"])
hos_interval_0_down_peak <- sum(hos_interval_0[which(hos_interval_0$date==hos_peakday_mean_interval_0),"down"])

p1 <- ggplot()+
  geom_bar(data=hos_interval_3,aes(x = date, y = mean/scale,fill=age,group=age), stat = "identity")+
  geom_errorbar(data=hos_ci_interval_3,aes(x = date, ymax = up/scale, ymin = down/scale),color=cl1,size=0.25)+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  geom_segment(aes(x=0,xend=xmax,
                   y=capacity_hos/scale,
                   yend=capacity_hos/scale),
               color="#b2182b", linetype="dashed",size=0.8)+
  geom_segment(aes(x=0,xend=xmax,
                   y=hos_interval_4_mean_peak/scale,
                   yend=hos_interval_4_mean_peak/scale),
               color="grey", linetype="dashed",size=0.8)+
  geom_ribbon(aes(x = c(0:xmax), 
                  ymin=hos_interval_4_down_peak/scale,
                  ymax=hos_interval_4_up_peak/scale),
              fill=cl2,alpha=0.5)+
  labs(x="Days since arrival of COVID-19 in China",y=paste("Hospital beds required (100,000)",sep=""))+
  scale_y_continuous(breaks=scales::pretty_breaks(5),labels = function(x) str_pad(x, 4, 'left', ' '),
                     sec.axis = sec_axis(~.*scale/capacity_hos*100,name = "Hospital occupancy (%)",
                                         breaks=scales::pretty_breaks(2),
                                         labels = function(x) str_pad(x, 5, 'right', ' ')))+
  mytheme+theme(panel.border = theme_border(type = c("right","bottom","left")))
p1

p2 <- ggplot()+
  geom_segment(aes(x=0,xend=xmax,
                   y=hos_interval_0_mean_peak/scale,
                   yend=hos_interval_0_mean_peak/scale),
               color="grey", linetype="solid",size=0.8)+
  geom_ribbon(aes(x = c(0:xmax),
                  ymin=hos_interval_0_down_peak/scale,
                  ymax=hos_interval_0_up_peak/scale),
              fill=cl2,alpha=0.5)+
  # geom_hline(yintercept=sum(hos_interval_0[which(hos_interval_0$date==peakday_mean_interval_0),"mean"])/scale,
  #            color="grey", linetype="dashed",size=0.8)+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(limits=c(430,520),breaks=scales::pretty_breaks(2),
                     sec.axis = sec_axis(~.*scale/capacity_hos*100,
                                         name = "Hospital occupancy (%)",
                                         breaks=c(1500,1600)))+
  labs(x=NULL,y=paste("Hospital",sep=""))+
  mytheme+theme(panel.border = theme_border(type = c("top","left","right")))


p_hos <- ggarrange(p2,p1,nrow=2,heights = c(1,3))
p_hos

# icu
xmin <- 0
xmax <- max(icu_interval_3$date)
icu_interval_4_mean_peak <- sum(icu_interval_4[which(icu_interval_4$date==icu_peakday_mean_interval_4),"mean"])
icu_interval_4_up_peak <- sum(icu_interval_4[which(icu_interval_4$date==icu_peakday_mean_interval_4),"up"])
icu_interval_4_down_peak <- sum(icu_interval_4[which(icu_interval_4$date==icu_peakday_mean_interval_4),"down"])
icu_interval_0_mean_peak <- sum(icu_interval_0[which(icu_interval_0$date==icu_peakday_mean_interval_0),"mean"])
icu_interval_0_up_peak <- sum(icu_interval_0[which(icu_interval_0$date==icu_peakday_mean_interval_0),"up"])
icu_interval_0_down_peak <- sum(icu_interval_0[which(icu_interval_0$date==icu_peakday_mean_interval_0),"down"])

p3 <- ggplot()+
  geom_bar(data=icu_interval_3,aes(x = date, y = mean/scale,fill=age,group=age), stat = "identity")+
  geom_errorbar(data=icu_ci_interval_3,aes(x = date, ymax = up/scale, ymin = down/scale),color=cl1,size=0.25)+
  geom_segment(aes(x=0,xend=xmax,
                   y=capacity_icu/scale,
                   yend=capacity_icu/scale),
               color="#b2182b", linetype="dashed",size=0.8)+
  geom_segment(aes(x=xmin,xend=xmax,
                   y=icu_interval_4_mean_peak/scale,
                   yend=icu_interval_4_mean_peak/scale),
                   color="grey", linetype="dashed",size=0.8)+
  geom_ribbon(aes(x = c(0:xmax), 
                  ymin=icu_interval_4_down_peak/scale,
                  ymax=icu_interval_4_up_peak/scale),
              fill=cl2,alpha=0.5)+
  scale_fill_manual(values = rev(c("#433983", "#20908c", "#34b778","#fde625")))+
  labs(x="Days since arrival of COVID-19 in China",y=paste("ICU beds required (100,000)",sep=""))+
  scale_y_continuous(limits=c(0,1.7),breaks=scales::pretty_breaks(9),labels = function(x) str_pad(x, 0, 'left', ' '),
                     sec.axis = sec_axis(~.*scale/capacity_icu*100,name = "ICU occupancy (%)",
                                         breaks=scales::pretty_breaks(3),
                                         labels = function(x) str_pad(x, 5, 'right', ' ')))+
  mytheme+theme(panel.border = theme_border(type = c("right","bottom","left")))
p3

p4 <- ggplot()+
  geom_segment(aes(x=xmin,xend=xmax,
                   y=icu_interval_0_mean_peak/scale,
                   yend=icu_interval_0_mean_peak/scale),
               color="grey", linetype="solid",size=0.8)+
  geom_ribbon(aes(x = c(0:xmax), 
                  ymin=icu_interval_0_down_peak/scale,
                  ymax=icu_interval_0_up_peak/scale),
              fill=cl2,alpha=0.5)+
  # geom_hline(yintercept=sum(icu_interval_0[which(icu_interval_0$date==icu_peakday_mean_interval_0),"mean"])/scale,
  #            color="grey", linetype="dashed",size=0.8)+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(limits=c(40.5,45.5),breaks=scales::pretty_breaks(3),labels = function(x) str_pad(x, 3, 'left', ' '),
                     sec.axis = sec_axis(~.*scale/capacity_icu*100,
                                         name = "icupital occupancy (%)",
                                         breaks=scales::pretty_breaks(2)))+
  labs(x=NULL,y=paste("Hospital",sep=""))+
  mytheme+theme(panel.border = theme_border(type = c("top","left","right")))
p4

p_icu <- ggarrange(p4,p3,nrow=2,heights = c(1,3))
p_icu
  
p <- ggarrange(p_hos,p_icu,ncol=2)
p

pdf(paste(filepath,"\\output\\fig3_a.pdf",sep=""),width=10,height=3.5)
print(p)
dev.off()

###hos
hos_interval_3_mean_peak <- sum(hos_interval_3[which(hos_interval_3$date==hos_peakday_mean_interval_3),"mean"])
hos_interval_3_up_peak <- sum(hos_interval_3[which(hos_interval_3$date==hos_peakday_mean_interval_3),"up"])
hos_interval_3_down_peak <- sum(hos_interval_3[which(hos_interval_3$date==hos_peakday_mean_interval_3),"down"])
round(hos_interval_3_mean_peak/capacity_hos*100,2)
round(hos_interval_3_up_peak/capacity_hos*100,2)
round(hos_interval_3_down_peak/capacity_hos*100,2)

round(hos_interval_0_mean_peak/capacity_hos,2)
round(hos_interval_0_up_peak/capacity_hos,2)
round(hos_interval_0_down_peak/capacity_hos,2)

###icu
icu_interval_3_mean_peak <- sum(icu_interval_3[which(icu_interval_3$date==icu_peakday_mean_interval_3),"mean"])
icu_interval_3_up_peak <- sum(icu_interval_3[which(icu_interval_3$date==icu_peakday_mean_interval_3),"up"])
icu_interval_3_down_peak <- sum(icu_interval_3[which(icu_interval_3$date==icu_peakday_mean_interval_3),"down"])
round(icu_interval_3_mean_peak/capacity_icu*100,2)
round(icu_interval_3_up_peak/capacity_icu*100,2)
round(icu_interval_3_down_peak/capacity_icu*100,2)

round(icu_interval_0_mean_peak/capacity_icu,2)
round(icu_interval_0_up_peak/capacity_icu,2)
round(icu_interval_0_down_peak/capacity_icu,2)
