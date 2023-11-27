rm(list = ls())

#### Packages ####
library(KRMr)
library(imager)
library(ggplot2)
library(tidyverse)
library(icesTAF)
library(doParallel)
library(TSdist)
library(dtw)


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

metrics.fast_slow.plot$lengthCat <- cut(metrics.fast_slow.plot$length, breaks=seq(from=20,to=50,by=5)*1e-2, right = FALSE)
metrics.model.plot$lengthCat <- cut(metrics.model.plot$length, breaks=seq(from=20,to=50,by=5)*1e-2, right = FALSE)

fast_slow.nspectra <- metrics.fast_slow.plot %>% 
                      group_by(freqChan,angleMatch,fish_id,species) %>% 
                      summarize(n.fast=mean(n.fast,na.rm=TRUE),
                                n.slow=mean(n.slow,na.rm=TRUE),
                                n.tot=(mean(n.fast,na.rm=TRUE)+mean(n.slow,na.rm=TRUE))/2)

metrics.model.plot$angleCat <- NA

idx.filt1 <- metrics.model.plot$angleMatch >= 245 & metrics.model.plot$angleMatch <= 295
metrics.model.plot$angleCat[idx.filt1] <- '245-295'
metrics.model.plot$angleCat[!idx.filt1] <- '65-115'

  
  


################################################################################
## plotting
################################################################################

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
  ylim(0,25)+
  ylab('Count (# per angle bin)')+
  xlab('Angle bin (1°)')

print(p)
dev.off()

# histogram of fish length
scaling_factor <- 1
png(file.path(myPath$figuresPath,'fish length histogram.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(summaryTab,ramping == 'slow'),aes(x=fish_length,fill=type))+
  geom_histogram()

print(p)
dev.off()

# consistency slow/fast ramping
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

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast min_species length.png'),
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=lengthCat,y=minDist,fill=species))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Minimum distance per angle bin (dB per kHz)')+
  xlab('Fish length category (cm)')+
  ylim(0,1)

print(p)
dev.off()

ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=lengthCat,y=meanTS.fast,fill=species))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Minimum distance per angle bin (dB per kHz)')+
  xlab('Fish length category (cm)')

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast TS.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanTS.slow))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Angle bin (1°)')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast TS_length species.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=lengthCat,y=meanTS.slow,fill=species))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Fish length category (cm)')+
  ylim(-60,-35)

print(p)
dev.off()

###############################################################################
# mean TS level betwee model and measurements
###############################################################################

windows()
ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=lengthCat,y=meanTS.mes,fill=species))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Mean TS (dB per kHz)')+
  xlab('Fish length category (cm)')+
  ylim(-70,-25)

windows()
ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=lengthCat,y=meanTS.model,fill=species))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~freqChan,ncol=1)+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Fish length category (cm)')+
  ylim(-70,-25)

scaling_factor <- 2
png(file.path(myPath$figuresPath,'TS_model vs measurements.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p1 <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanTS.model))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(freqChan~angleCat,scales='free_x')+
  ylab(expression('TS (dB re 1m'^'2'*')'))+
  xlab('Angle bin (1°)')+
  ylim(-60,-25)+
  scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                   labels = c("65",'90','115','245','270','295'))+
  ggtitle('Model')

p2 <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanTS.mes))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(freqChan~angleCat,scales='free_x')+
  ylab('Mean TS (dB per kHz)')+
  xlab('Angle bin (1°)')+
  ylim(-60,-25)+
  scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                   labels = c("65",'90','115','245','270','295'))+
  ggtitle('Measurements')

p.all <- grid.arrange(p1, p2, nrow = 1)

print(p.all)
dev.off()
###############################################################################
# model vs measurements
###############################################################################



###############################################################################
# model vs measurements
###############################################################################

scaling_factor <- 2
png(file.path(myPath$figuresPath,'model vs measurements minDist.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=minDist,fill=ramping))+
  geom_boxplot()+
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
      scale_x_discrete(breaks= c("65",'90','115','245','270','295'),#unique(SA.station$time_day)[1:10],
                       labels = c("65",'90','115','245','270','295'))+
      ylim(0,2)

print(p)
dev.off()










scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast meanDistTrunc.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanDistTrunc25))+
  theme_bw()+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Mean distance per angle bin (dB per kHz)')+
  xlab('Angle bin (1°)')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,1.5)

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'slow fast sd.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.fast_slow.plot,freqChan != 'all'),aes(x=as.factor(angleMatch),y=meanDistTrunc25))+
  theme_bw()+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~freqChan,ncol=1)+
  ylab('Standard deviation per angle bin (dB per kHz)')+
  xlab('Angle bin (1°)')+
  scale_x_discrete(breaks= c("0",'90','180','270','360'),#unique(SA.station$time_day)[1:10],
                   labels = c("0",'90','180','270','360'))+
  ylim(0,0.5)

print(p)
dev.off()

A <- subset(metrics.fast_slow.plot,freqChan == '200khz' & angleMatch >= 89 & angleMatch <= 89)# & species =='saithe'

windows()
ggplot()+
  geom_point(data=A$meanTS.slow,
             aes(x=(bodyV),y=meanDistAngleBin.slow,col='slow'))+
  geom_point(data=A,
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

# matching between model and data
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

scaling_factor <- 1.5
png(file.path(myPath$figuresPath,'model vs fast mean_fish length_90deg.png'), 
    width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.model.plot,freqChan != 'all' & angleMatch == 90),aes(x=fish_length,y=mean,col=ramping))+
  geom_point()+
  facet_wrap(~freqChan)+
  ggtitle('90 degree angle')

print(p)
dev.off()

