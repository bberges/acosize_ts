library(ggplot2)
library(tidyverse)
library(ggbiplot)
library(mgcv)
library(ale)

rm(list = ls())

path <- 'G:/acosize_ts'

try(setwd(path),silent=TRUE)

species <- 'SAI'

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     rPath=file.path(".","R"))


predTab <- read.csv(file = file.path(myPath$resultsPath,paste0('pred_',species,'.csv')))
predTab$angleGroup <- cut(predTab$angle,
                         seq(min(predTab$angle),360,by=36),
                         include.lowest = T)
predTab$angleGroup <- as.factor(predTab$angleGroup)

p <- ggplot(predTab,aes(x=length,y=length_pred,group=angleGroup))+
      stat_summary(fun.y = mean,
                   fun.min = function(x) mean(x) - sd(x),
                   fun.max = function(x) mean(x) + sd(x),
                   geom = "pointrange")+
      ylab('predictions')+
      xlab('targets')+
      ggtitle(species)+
      geom_smooth(method='lm', formula= y~x)+
      facet_wrap(~angleGroup)

ggsave(file.path(myPath$figuresPath,paste0('pred_singleTarget_',species,'.png')),
       p,
       width = 170,
       height = 150,
       units = c("mm"),
       dpi = 300)


predAngleTab <- read.csv(file = file.path(myPath$resultsPath,paste0('pred_',species,'_angle.csv')))
predAngleTab$length <- as.numeric(predAngleTab$length)
predAngleTab$length_pred <- as.numeric(predAngleTab$length_pred)
predAngleTab <- subset(predAngleTab,!is.na(length) | !is.na(length_pred))


windows()
ggplot(predAngleTab,aes(x=as.factor(length),y=length_pred))+
  geom_boxplot()


windows()
ggplot(predAngleTab,aes(x=length,y=length_pred,group=length))+
  stat_summary(fun.y = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "pointrange")+
  ylab('Angle predictions')+
  xlab('Angle targets')


windows()
ggplot(predAngleTab,aes(x=length*180/pi,y=length_pred*180/pi,group=length*180/pi))+
  stat_summary(fun.y = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "pointrange")+
  ylab('Angle predictions')+
  xlab('Angle targets')
