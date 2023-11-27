rm(list = ls())

#### Packages ####
library(KRMr)
library(imager)
library(ggplot2)
library(tidyverse)
library(icesTAF)
library(doParallel)
library(gridExtra)


#path <- 'J:/git/heras_index_kbwot'
path <- 'C:/git/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     rPath=file.path(".","R"))

distanceMeasure <- 'eucl'
windowing <- 'variableWindow' # variableWindow fixedWindow
dir.create(file.path(myPath$figuresPath,paste0(windowing,'_',distanceMeasure)),showWarnings = FALSE)
myPath$figuresPath <- file.path(myPath$figuresPath,paste0(windowing,'_',distanceMeasure))

################################################################################
## load objects
################################################################################

summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]
summaryTab <- summaryTab[summaryTab$orientation == 'ventral',]

load(file.path(myPath$resultsPath,'fish_dimensions.RData'))

load(file.path(myPath$resultsPath,paste0('acosize_TS_summary results_',
                                         windowing, '_',
                                         distanceMeasure,'.RData')))

uniqueAnimals <- unique(summaryTab$fish_id)

################################################################################
## formatting tables
################################################################################
metrics.model.all     <- metrics.model.all[!is.na(metrics.model.all$freqChan),]
metrics.slow_fast.all <- metrics.slow_fast.all[!is.na(metrics.slow_fast.all$freqChan),]

tabJoin <- subset(summaryTab,ramping == 'fast') %>% select(-c('ramping'))
metrics.model.plot <- left_join(metrics.model.all,metrics.shape,by='fish_id')

metrics.fast_slow.plot <- left_join(metrics.slow_fast.all,metrics.shape,by='fish_id')

fast_slow.nspectra <- metrics.fast_slow.plot %>% 
                      group_by(freqChan,angleMatch,fish_id,species) %>% 
                      summarize(fast=mean(n.fast,na.rm=TRUE),
                                slow=mean(n.slow,na.rm=TRUE)) %>%
                      pivot_longer(!freqChan & !angleMatch & !fish_id & !species,names_to = 'ramping',values_to = 'N')

fast_slow.nspectra.red <- subset(fast_slow.nspectra,freqChan != 'all') %>%
                          group_by(freqChan,angleMatch,species) %>%
                          summarize(N=mean(N,na.rm=TRUE))

metrics.model.plot$angleCat <- NA

idx.filt1 <- metrics.model.plot$angleMatch >= 245 & metrics.model.plot$angleMatch <= 295
metrics.model.plot$angleCat[idx.filt1] <- '245-295'
metrics.model.plot$angleCat[!idx.filt1] <- '65-115'
metrics.model.plot$angleCat <- factor(metrics.model.plot$angleCat,levels=c('65-115','245-295'))

metrics.fast_slow.plot$lengthCat <- cut(metrics.fast_slow.plot$length, breaks=seq(from=20,to=50,by=2)*1e-2, right = FALSE)
metrics.model.plot$lengthCat <- cut(metrics.model.plot$length, breaks=seq(from=20,to=50,by=2)*1e-2, right = FALSE)

################################################################################
## plotting
################################################################################

############################
# WGFSAT plots
############################

myPath$figuresPath <- file.path(myPath$figuresPath,'WGFAST')

scaling_factor <- 2
png(file.path(myPath$figuresPath,'TS_model vs measurements_70khz.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p1 <- ggplot(subset(metrics.model.plot,freqChan != 'all' & freqChan == '70khz'),aes(x=as.factor(angleMatch),y=meanTS.model,fill=species))+
      theme_bw()+
      geom_boxplot(fatten = NULL,lwd=0.1)+
      geom_boxplot(aes(color=species))+
      facet_wrap(~angleCat,scales='free_x')+
      ylab(expression('TS (dB re 1m'^'2'*')'))+
      xlab('Angle bin (1°)')+
      ylim(-70,-25)+
      scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                       labels = c("65",'90','115','245','270','295'))+
      theme(axis.text.x = element_text(angle = 90))+
      ggtitle('Model - 70 kHz')

p2 <- ggplot(subset(metrics.model.plot,freqChan != 'all' & freqChan == '70khz'),aes(x=as.factor(angleMatch),y=meanTS.mes,fill=species))+
      theme_bw()+
      geom_boxplot()+
      geom_boxplot(aes(color=species))+
      facet_wrap(~angleCat,scales='free_x')+
      ylab(expression('TS (dB re 1m'^'2'*')'))+
      xlab('Angle bin (1°)')+
      ylim(-70,-25)+
      scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                       labels = c("65",'90','115','245','270','295'))+
      theme(axis.text.x = element_text(angle = 90))+
      ggtitle('Measurements - 70 kHz')

p.70khz <- grid.arrange(p1, p2, ncol = 1)

print(p.70khz)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'TS_model vs measurements_200khz.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p1 <- ggplot(subset(metrics.model.plot,freqChan != 'all' & freqChan == '200khz'),aes(x=as.factor(angleMatch),y=meanTS.model,fill=species))+
  theme_bw()+
  geom_boxplot(fatten = NULL,lwd=0.1)+
  geom_boxplot(aes(color=species))+
  facet_wrap(~angleCat,scales='free_x')+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Angle bin (1°)')+
  ylim(-70,-25)+
  scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                   labels = c("65",'90','115','245','270','295'))+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Model - 200 kHz')

p2 <- ggplot(subset(metrics.model.plot,freqChan != 'all' & freqChan == '200khz'),aes(x=as.factor(angleMatch),y=meanTS.mes,fill=species))+
  theme_bw()+
  geom_boxplot()+
  geom_boxplot(aes(color=species))+
  facet_wrap(~angleCat,scales='free_x')+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Angle bin (1°)')+
  ylim(-70,-25)+
  scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                   labels = c("65",'90','115','245','270','295'))+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Measurements - 200 kHz')

p.200khz <- grid.arrange(p1, p2, ncol = 1)

print(p.200khz)
dev.off()

# N spectra
scaling_factor <- 2
png(file.path(myPath$figuresPath,'n spectra angle bins.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(data=subset(fast_slow.nspectra,freqChan != 'all'),aes(x=as.factor(angleMatch),y=N,fill=freqChan))+
      geom_boxplot(fatten = NULL,outlier.shape = NA)+
      geom_boxplot(aes(col=freqChan),outlier.shape = NA)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~species,ncol=1)+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,25)+
  ylab('Count (# per angle bin)')+
  xlab('Angle bin (1 degree)')

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'n spectra angle bins_2.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(fast_slow.nspectra.red,aes(x=angleMatch,y=N,col=freqChan))+
      geom_point()+
      facet_wrap(~species,ncol=1)+
      ylim(0,25)+
      ylab('Count (# per angle bin)')+
      xlab('Angle bin (1 degree)')

print(p)
dev.off()




scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=minDist))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,1)

print(p)
dev.off()


scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs measurements minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=minDist,fill=ramping))+
  geom_boxplot(fatten = NULL,outlier.shape = NA)+
  geom_boxplot(aes(col=ramping),outlier.shape = NA)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_grid(freqChan~angleCat,scales='free_x')+
  ylab('Minimum distance per angle bin (dB per kHz)')+
  xlab('Angle bin (1°)')+
  scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                   labels = c("65",'90','115','245','270','295'))+
  ylim(0,2)

print(p)
dev.off()


scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs measurements minDist_length species.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(lengthCat),y=minDist,fill=species))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_grid(freqChan~angleCat,scales='free_x')+
  ylab('Minimum distance per angle bin (dB per kHz)')+
  xlab('Angle bin (1°)')+
  ylim(0,2)

print(p)
dev.off()


scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs measurements TS(L) 89-91.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

temp.plot <- subset(metrics.model.plot,freqChan != 'all' & angleMatch >= 89 & angleMatch <= 91)
p1 <- ggplot(temp.plot,aes(x=length,y=meanTS.mes,col=species))+#,col=ramping
  geom_point()+
  geom_errorbar(aes(x=length,ymin=meanTS.mes-sdTS.mes, ymax=meanTS.mes+sdTS.mes))+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  ylim(-60,-30)+
  ggtitle('Measurements - TS vs length 89-91 degree')


p2 <- ggplot(temp.plot,aes(x=length,y=meanTS.model,col=species))+#,col=ramping
  geom_point()+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  ylim(-60,-30)+
  ggtitle('Model - TS vs length 89-91 degree')


print(grid.arrange(p2, p1, ncol = 1))
dev.off()


scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs measurements TS(L) 63-65.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

temp.plot <- subset(metrics.model.plot,freqChan != 'all' & angleMatch >= 63 & angleMatch <= 65)
p1 <- ggplot(temp.plot,aes(x=length,y=meanTS.mes,col=species))+#,col=ramping
      geom_point()+
      geom_errorbar(aes(x=length,ymin=meanTS.mes-sdTS.mes, ymax=meanTS.mes+sdTS.mes))+
      facet_wrap(~freqChan)+
      stat_smooth(method = "lm")+
      ylab(expression('TS (dB re 1m'^'2'*')'))+
      ylim(-70,-30)+
      ggtitle('Measurements - TS vs length 63-65 degree')


p2 <- ggplot(temp.plot,aes(x=length,y=meanTS.model,col=species))+#,col=ramping
  geom_point()+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  ylim(-70,-30)+
  ggtitle('Model - TS vs length 63-65 degree')


print(grid.arrange(p2, p1, ncol = 1))
dev.off()





temp.plot <- subset(metrics.model.plot,freqChan != 'all' & angleMatch >= 89 & angleMatch <= 91)
p1 <- ggplot(temp.plot,aes(x=length,y=meanTS.mes,col=species))+#,col=ramping
      geom_point()+
      geom_errorbar(aes(x=length,ymin=meanTS.mes-sdTS.mes, ymax=meanTS.mes+sdTS.mes))+
      facet_wrap(~freqChan)+
      stat_smooth(method = "lm")+
      ylab(expression('TS (dB re 1m'^'2'*')'))+
      ylim(-70,-30)+
      ggtitle('TS vs length - measurements')


p2 <- ggplot(temp.plot,aes(x=length,y=meanTS.model,col=species))+#,col=ramping
      geom_point()+
      facet_wrap(~freqChan)+
      stat_smooth(method = "lm")+
      ylab(expression('TS (dB re 1m'^'2'*')'))+
      ylim(-70,-30)+
      ggtitle('TS vs length - model')

windows()
grid.arrange(p1, p2, ncol = 1)


temp.plot <- subset(metrics.model.plot,freqChan != 'all' & angleMatch >= 89 & angleMatch <= 91)
p1 <- ggplot(temp.plot,aes(x=length,y=meanTS.mes,col=species))+#,col=ramping
  geom_point()+
  geom_errorbar(aes(x=length,ymin=meanTS.mes-sdTS.mes, ymax=meanTS.mes+sdTS.mes))+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  ylim(-70,-30)+
  ggtitle('TS vs length - measurements')


p2 <- ggplot(temp.plot,aes(x=length,y=meanTS.model,col=species))+#,col=ramping
  geom_point()+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  ylim(-60,-30)+
  ggtitle('TS vs length - model')

windows()
grid.arrange(p1, p2, ncol = 1)


############################
# number of spectra
############################
scaling_factor <- 2
png(file.path(myPath$figuresPath,'n spectra angle bins_fast.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(data=subset(fast_slow.nspectra,!is.na(n.fast) & freqChan != '70khz'),aes(x=as.factor(angleMatch),y=n.fast))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~species,ncol=1)+
  ggtitle('n spectra fast ramping')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,30)+
  ylab('Count (# per angle bin)')+
  xlab('Angle bin (1 degree)')

print(p)
dev.off()



scaling_factor <- 2
png(file.path(myPath$figuresPath,'n spectra angle bins_slow.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(data=subset(fast_slow.nspectra,!is.na(n.fast) & freqChan != '70khz'),aes(x=as.factor(angleMatch),y=n.slow))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~species,ncol=1)+
  ggtitle('n spectra slow ramping')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,30)+
  ylab('Count (# per angle bin)')+
  xlab('Angle bin (1 degree)')

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'n spectra angle bins_tot.png'),
    width = 10*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(data=subset(fast_slow.nspectra,!is.na(n.fast) & freqChan != '70khz'),aes(x=as.factor(angleMatch),y=n.tot))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~species,ncol=1)+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,30)+
  ylab('Count (# per angle bin)')+
  xlab('Angle bin (1 degree)')

print(p)
dev.off()

############################
# histogram of fish length
############################
scaling_factor <- 1
png(file.path(myPath$figuresPath,'fish length histogram.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(summaryTab,ramping == 'slow'),aes(x=fish_length,fill=type))+
  geom_histogram()

print(p)
dev.off()

############################
# TS levels
############################
temp.plot <- subset(metrics.fast_slow.plot,freqChan == '200khz' & angleMatch >= 89 & angleMatch <= 89)# & species =='saithe'

windows()
ggplot()+
  geom_point(data=temp.plot,
             aes(x=(bodyV),y=meanDistAngleBin.slow,col='slow'))+
  geom_point(data=temp.plot,
             aes(x=(bodyV),y=meanDistAngleBin.fast,col='fast'))+
  facet_grid(~species,scales='free')

windows()
ggplot(subset(metrics.fast_slow.plot,freqChan != 'all' & angleMatch >= 88 & angleMatch <= 92),
       aes(x=as.factor(length),y=meanTS.fast))+
  geom_point()+
  facet_grid(freqChan~species)

windows()
ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanTS.fast))+
  geom_boxplot()+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))

temp.plot <- subset(metrics.model.plot,freqChan != 'all' & angleMatch >= 89 & angleMatch <= 91)
windows()
ggplot(temp.plot,aes(x=length,y=meanTS.mes,col=species))+#,col=ramping
  geom_point()+
  geom_errorbar(aes(x=length,ymin=meanTS.mes-sdTS.mes, ymax=meanTS.mes+sdTS.mes))+
  facet_wrap(~freqChan)+
  stat_smooth(method = "lm")+
  ggtitle('TS vs length')

############################
# consistency slow/fast ramping
############################
scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanDistTrunc25))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast sd.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=sdDist))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,0.5)

print(p)
dev.off()


scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast mean_fish length.png'),
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(length),y=meanDist,fill=species))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast min_fish length.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(fish_length),y=minDist,fill=type))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

###################################
# model match
###################################
scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast meanDistTrunc25.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanDistTrunc25,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast mean.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=mean,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=minDist,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast sd.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=sdDist,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,0.5)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast meanDistTrunc25_fish length.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(fish_length),y=meanDistTrunc25,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,2)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs fast sd_fish length.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(fish_length),y=sdDist,fill=ramping))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylim(0,0.5)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=minDist))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Minimum distance per angle bin (dB per kHz)')+
  xlab('Angle bin (1°)')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,1)

print(p)
dev.off()

