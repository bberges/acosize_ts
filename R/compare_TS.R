rm(list = ls())

#### Packages ####
library(KRMr)
library(imager)
library(ggplot2)
library(tidyverse)
library(icesTAF)
library(doParallel)
library(TSdist)


#path <- 'J:/git/heras_index_kbwot'
path <- 'G:/git/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"))

sourceDir(myPath$mainPath)

summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]
summaryTab <- summaryTab[summaryTab$orientation == 'ventral',]

uniqueAnimals <- unique(summaryTab$fish_id)

flagFirstAnimal <- TRUE
for(idxAnimal in uniqueAnimals){
  print(idxAnimal)
  KRMCurrent <- read.csv(file.path(myPath$resultsPath,'KRM',
                                   paste0('krm_',idxAnimal,'.csv')))
  
  currentTab <- summaryTab[summaryTab$fish_id == idxAnimal,]

  flagFirst <- TRUE
  for(idxRamping in unique(currentTab$ramping)){
    measuName <- paste0(unique(currentTab$fish_id),'_',
                        unique(currentTab$mode),'_',
                        idxRamping,'_',
                        unique(currentTab$orientation),'_',
                        unique(currentTab$resolution))
    
    ###################################
    # load 200khz results
    ###################################
    measuName.200khz  <- paste0(measuName,'_200khz')
    TS.cor.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                        paste0(measuName.200khz,'_TSCorr.csv')),check.names=FALSE)
    TS.cor.200khz$angleMatch <- TS.cor.200khz$angle-270
    TS.cor.200khz$angleMatch[TS.cor.200khz$angleMatch < 0] <- TS.cor.200khz$angleMatch[TS.cor.200khz$angleMatch < 0]+360
    TS.cor.200khz$freqChan <- '200khz'
    TS.cor.200khz$ramping <- idxRamping
    
    TS.kde.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                        paste0(measuName.200khz,'_TSkde.csv')),check.names=FALSE)
    TS.kde.200khz <- TS.kde.200khz %>% pivot_longer(!N_spectra & !angle & !TSmesh,names_to = "frequency", values_to = "density")
    TS.kde.200khz$frequency <- as.numeric(TS.kde.200khz$frequency)
    TS.kde.200khz$angleMatch <- TS.kde.200khz$angle-270
    TS.kde.200khz$angleMatch[TS.kde.200khz$angleMatch < 0] <- TS.kde.200khz$angleMatch[TS.kde.200khz$angleMatch < 0]+360
    TS.kde.200khz$freqChan <- '200khz'
    TS.kde.200khz$ramping <- idxRamping
    
    TS.stats.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                          paste0(measuName.200khz,'_TSStats.csv')),check.names=FALSE)
    TS.stats.200khz <- TS.stats.200khz %>% pivot_longer(!N_spectra & !angle & !percentiles,names_to = "frequency", values_to = "TS")
    TS.stats.200khz$frequency <- as.numeric(TS.stats.200khz$frequency)
    TS.stats.200khz$angleMatch <- TS.stats.200khz$angle-270
    TS.stats.200khz$angleMatch[TS.stats.200khz$angleMatch < 0] <- TS.stats.200khz$angleMatch[TS.stats.200khz$angleMatch < 0]+360
    TS.stats.200khz$freqChan <- '200khz'
    TS.stats.200khz$ramping <- idxRamping
    
    ###################################
    # load 70khz results
    ###################################
    measuName.70khz   <- paste0(measuName,'_70khz')
    TS.cor.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                       paste0(measuName.70khz,'_TSCorr.csv')),check.names=FALSE)
    TS.cor.70khz$angleMatch <- TS.cor.70khz$angle-270
    TS.cor.70khz$angleMatch[TS.cor.70khz$angleMatch < 0] <- TS.cor.70khz$angleMatch[TS.cor.70khz$angleMatch < 0]+360
    TS.cor.70khz$freqChan <- '70khz'
    TS.cor.70khz$ramping <- idxRamping
    
    TS.kde.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                       paste0(measuName.70khz,'_TSkde.csv')),check.names=FALSE)
    TS.kde.70khz <- TS.kde.70khz %>% pivot_longer(!N_spectra & !angle & !TSmesh,names_to = "frequency", values_to = "density")
    TS.kde.70khz$frequency <- as.numeric(TS.kde.70khz$frequency)
    TS.kde.70khz$angleMatch <- TS.kde.70khz$angle-270
    TS.kde.70khz$angleMatch[TS.kde.70khz$angleMatch < 0] <- TS.kde.70khz$angleMatch[TS.kde.70khz$angleMatch < 0]+360
    TS.kde.70khz$freqChan <- '70khz'
    TS.kde.70khz$ramping <- idxRamping
    
    TS.stats.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                         paste0(measuName.70khz,'_TSStats.csv')),check.names=FALSE)
    TS.stats.70khz <- TS.stats.70khz %>% pivot_longer(!N_spectra & !angle & !percentiles,names_to = "frequency", values_to = "TS")
    TS.stats.70khz$frequency <- as.numeric(TS.stats.70khz$frequency)
    TS.stats.70khz$angleMatch <- TS.stats.70khz$angle-270
    TS.stats.70khz$angleMatch[TS.stats.70khz$angleMatch < 0] <- TS.stats.70khz$angleMatch[TS.stats.70khz$angleMatch < 0]+360
    TS.stats.70khz$freqChan <- '70khz'
    TS.stats.70khz$ramping <- idxRamping
    
    ###################################
    # combine tables
    ###################################
    if(flagFirst){
      TS.cor.all <- rbind(TS.cor.70khz,TS.cor.200khz)
      TS.stats.all <- rbind(TS.stats.70khz,TS.stats.200khz)
      TS.kde.all <- rbind(TS.kde.70khz,TS.kde.200khz)
      flagFirst <- FALSE
    }else{
      TS.cor.all <- rbind(TS.cor.all,TS.cor.70khz,TS.cor.200khz)
      TS.stats.all <- rbind(TS.stats.all,TS.stats.70khz,TS.stats.200khz)
      TS.kde.all <- rbind(TS.kde.all,TS.kde.70khz,TS.kde.200khz)
    }
  }
  
  KRMCurrent.filt <- KRMCurrent[KRMCurrent$frequency %in% unique(TS.stats.all$frequency),]
  uniqueAngles.model <- unique(KRMCurrent.filt$theta)
  
  freq.70khz <- unique(TS.stats.all$frequency[TS.stats.all$freqChan == '70khz'])
  freq.200khz <- unique(TS.stats.all$frequency[TS.stats.all$freqChan == '200khz'])
  
  idx70khz <- unique(TS.stats.all$frequency) %in% freq.70khz
  idx200khz <- unique(TS.stats.all$frequency) %in% freq.200khz

  TS.stats.wide <- TS.stats.all %>% group_by(angle,percentiles,frequency,angleMatch,ramping) %>% summarize(TS=TS) %>% pivot_wider(names_from=frequency,values_from=TS)
  TS.stats.wide <- TS.stats.all %>% group_by(angle,percentiles,frequency,angleMatch,ramping) %>% summarize(TS=TS) %>% pivot_wider(names_from=frequency,values_from=TS)
  uniqueAngles <- unique(TS.stats.wide$angleMatch)
  
  
  diff.fast_slow <- as.data.frame(matrix(ncol = length(unique(TS.stats.all$frequency))+1, nrow = length(uniqueAngles)))
  colnames(diff.fast_slow) <- c('angleMatch',unique(TS.stats.all$frequency))
  
  metrics.fast_slow <- as.data.frame(matrix(ncol = 8, nrow = length(uniqueAngles)*3))
  colnames(metrics.fast_slow) <- c('angleMatch','freqChan','quant90','quant50','quant10','mean','cor','dist.cor')

  diff.model <- as.data.frame(matrix(ncol = length(unique(TS.stats.all$frequency))+2, nrow = length(uniqueAngles.model)*2))
  colnames(diff.model) <- c('angleMatch','ramping',unique(TS.stats.all$frequency))
  
  metrics.model <- as.data.frame(matrix(ncol = 9, nrow = length(uniqueAngles.model)*2*3))
  colnames(metrics.model) <- c('angleMatch','ramping','freqChan','quant90','quant50','quant10','mean','cor','dist.cor')

  counter.fast_slow <- 1
  counter.model <- 1
  for(idxAngle in uniqueAngles){
    TS.stats.wide.filter <- subset(TS.stats.wide,angleMatch == idxAngle & percentiles == 50)
    TS.stats.wide.filter <- TS.stats.wide.filter %>% ungroup() %>% select(-c('angle','percentiles','angleMatch'))
    
    #plot.df <- subset(TS.stats.wide.filter,ramping=='fast') %>% pivot_longer(!ramping,names_to = 'frequency',values_to = 'TS')
    #plot.df$frequency <- as.numeric(plot.df$frequency)
    #ggplot()+
    #  geom_line(data=KRMCurrent.current,aes(x=frequency,y=TS),col='blue')+
    #  geom_line(data=plot.df,aes(x=frequency,y=TS),col='red')
    
    ts.slow <- as.numeric(TS.stats.wide.filter[TS.stats.wide.filter$ramping == 'slow',] %>% select(-c('ramping')))
    ts.fast <- as.numeric(TS.stats.wide.filter[TS.stats.wide.filter$ramping == 'fast',] %>% select(-c('ramping')))
    myVec.fast_slow <- ts.fast-ts.slow
    
    diff.fast_slow[counter.fast_slow,2:(length(unique(TS.stats.all$frequency))+1)] <- myVec.fast_slow
    diff.fast_slow$angleMatch[counter.fast_slow] <- idxAngle
    
    metrics.fast_slow$angleMatch[counter.fast_slow] <- idxAngle
    metrics.fast_slow$quant90[counter.fast_slow] <- quantile(abs(myVec.fast_slow),c(0.9),na.rm=TRUE)
    metrics.fast_slow$quant50[counter.fast_slow] <- quantile(abs(myVec.fast_slow),c(0.5),na.rm=TRUE)
    metrics.fast_slow$quant10[counter.fast_slow] <- quantile(abs(myVec.fast_slow),c(0.1),na.rm=TRUE)
    metrics.fast_slow$mean[counter.fast_slow]    <- mean(abs(myVec.fast_slow),na.rm=TRUE)
    metrics.fast_slow$cor[counter.fast_slow]       <- cor(ts.slow,ts.fast)
    metrics.fast_slow$dist.cor[counter.fast_slow]  <- TSDistances(ts.slow, ts.fast, distance="cor")
    metrics.fast_slow$freqChan[counter.fast_slow] <- 'all'
    counter.fast_slow <- counter.fast_slow+1
    
    metrics.fast_slow$angleMatch[counter.fast_slow] <- idxAngle
    metrics.fast_slow$quant90[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx70khz]),c(0.9),na.rm=TRUE)
    metrics.fast_slow$quant50[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx70khz]),c(0.5),na.rm=TRUE)
    metrics.fast_slow$quant10[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx70khz]),c(0.1),na.rm=TRUE)
    metrics.fast_slow$mean[counter.fast_slow]    <- mean(abs(myVec.fast_slow[idx70khz]),na.rm=TRUE)
    metrics.fast_slow$cor[counter.fast_slow]       <- cor(ts.slow[idx70khz],ts.fast[idx70khz])
    metrics.fast_slow$dist.cor[counter.fast_slow]  <- TSDistances(ts.slow[idx70khz], ts.fast[idx70khz], distance="cor")
    metrics.fast_slow$freqChan[counter.fast_slow] <- '70khz'
    counter.fast_slow <- counter.fast_slow+1
    
    metrics.fast_slow$angleMatch[counter.fast_slow] <- idxAngle
    metrics.fast_slow$quant90[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx200khz]),c(0.9),na.rm=TRUE)
    metrics.fast_slow$quant50[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx200khz]),c(0.5),na.rm=TRUE)
    metrics.fast_slow$quant10[counter.fast_slow] <- quantile(abs(myVec.fast_slow[idx200khz]),c(0.1),na.rm=TRUE)
    metrics.fast_slow$mean[counter.fast_slow]    <- mean(abs(myVec.fast_slow[idx200khz]),na.rm=TRUE)
    metrics.fast_slow$cor[counter.fast_slow]       <- cor(ts.slow[idx200khz],ts.fast[idx200khz])
    metrics.fast_slow$dist.cor[counter.fast_slow]  <- TSDistances(ts.slow[idx200khz], ts.fast[idx200khz], distance="cor")
    metrics.fast_slow$freqChan[counter.fast_slow] <- '200khz'
    counter.fast_slow <- counter.fast_slow+1

    if(idxAngle %in% unique(KRMCurrent.filt$theta)){
      KRMCurrent.current <- subset(KRMCurrent.filt,theta == idxAngle) %>% dplyr::select(c('frequency','TS'))# %>% pivot_longer(names_from=frequency,values_from=TS)
      ts.model <- KRMCurrent.current$TS
      myVec.model.fast <- ts.fast-ts.model
      myVec.model.slow <- ts.slow-ts.model

      diff.model[counter.model,3:(length(unique(TS.stats.all$frequency))+2)]    <- myVec.model.fast
      diff.model[counter.model+1,3:(length(unique(TS.stats.all$frequency))+2)]  <- myVec.model.slow
      diff.model$angleMatch[counter.model:(counter.model+1)] <- idxAngle
      diff.model$ramping[counter.model] <- 'fast'
      diff.model$ramping[counter.model+1] <- 'slow'
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'fast'
      metrics.model$freqChan[counter.model] <- 'all'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.fast),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.fast),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.fast),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.fast),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.fast,ts.model)
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.fast, ts.model, distance="cor")
      counter.model <- counter.model+1
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'slow'
      metrics.model$freqChan[counter.model] <- 'all'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.slow),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.slow),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.slow),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.slow),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.slow,ts.model)
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.slow, ts.model, distance="cor")
      counter.model <- counter.model+1
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'fast'
      metrics.model$freqChan[counter.model] <- '70khz'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.fast[idx70khz]),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.fast[idx70khz]),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.fast[idx70khz]),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.fast[idx70khz]),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.fast[idx70khz],ts.model[idx70khz])
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.fast[idx70khz], ts.model[idx70khz], distance="cor")
      counter.model <- counter.model+1
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'slow'
      metrics.model$freqChan[counter.model] <- '70khz'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.slow[idx70khz]),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.slow[idx70khz]),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.slow[idx70khz]),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.slow[idx70khz]),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.slow[idx70khz],ts.model[idx70khz])
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.slow[idx70khz], ts.model[idx70khz], distance="cor")
      counter.model <- counter.model+1
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'fast'
      metrics.model$freqChan[counter.model] <- '200khz'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.fast[idx200khz]),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.fast[idx200khz]),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.fast[idx200khz]),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.fast[idx200khz]),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.fast[idx200khz],ts.model[idx200khz])
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.fast[idx200khz], ts.model[idx200khz], distance="cor")
      counter.model <- counter.model+1
      
      metrics.model$angleMatch[counter.model] <- idxAngle
      metrics.model$ramping[counter.model] <- 'slow'
      metrics.model$freqChan[counter.model] <- '200khz'
      metrics.model$quant90[counter.model]   <- quantile(abs(myVec.model.slow[idx200khz]),c(0.9),na.rm=TRUE)
      metrics.model$quant50[counter.model]   <- quantile(abs(myVec.model.slow[idx200khz]),c(0.5),na.rm=TRUE)
      metrics.model$quant10[counter.model]   <- quantile(abs(myVec.model.slow[idx200khz]),c(0.1),na.rm=TRUE)
      metrics.model$mean[counter.model]      <- mean(abs(myVec.model.slow[idx200khz]),na.rm=TRUE)
      metrics.model$cor[counter.model]       <- cor(ts.slow[idx200khz],ts.model[idx200khz])
      metrics.model$dist.cor[counter.model]  <- TSDistances(ts.slow[idx200khz], ts.model[idx200khz], distance="cor")
      counter.model <- counter.model+1
    }
  }
  
  metrics.model$fish_id <- idxAnimal
  diff.model$fish_id <- idxAnimal
  metrics.fast_slow$fish_id <- idxAnimal
  diff.fast_slow$fish_id <- idxAnimal
  
  if(flagFirstAnimal){
    metrics.model.all <- metrics.model
    diff.model.all <- diff.model
    metrics.fast_slow.all <- metrics.fast_slow
    diff.fast_slow.all <- diff.fast_slow
    flagFirstAnimal <- FALSE
  }else{
    metrics.model.all <- rbind(metrics.model.all,metrics.model)
    diff.model.all <- rbind(diff.model.all,diff.model)
    metrics.fast_slow.all <- rbind(metrics.fast_slow.all,metrics.fast_slow)
    diff.fast_slow.all <- rbind(diff.fast_slow.all,diff.fast_slow)
  }
  
  #windows()
  #ggplot(data=diff.fast_slow,aes(x=angleMatch))+
  #  geom_line(aes(y=quant50))+
  #  geom_ribbon(aes(ymin=quant10,ymax=quant90))+
  #  geom_line(aes(x=angleMatch,y=mean),col='red')
  
  #ggplot(diff.fast_slow,aes(x=as.factor(angleMatch),y=diff.stats))+
  #  geom_boxplot(outlier.shape = NA)+
  #  coord_cartesian(ylim = quantile(diff.fast_slow$diff.stats, c(0.1, 0.9),na.rm=TRUE))
}

save(metrics.model.all,
     diff.model.all,
     metrics.fast_slow.all,
     diff.fast_slow.all,
     file = file.path(myPath$resultsPath,'acosize_TS_summary results.RData'))

tabJoin <- subset(summaryTab,ramping == 'fast') %>% select(-c('ramping'))
metrics.model.plot <- left_join(metrics.model.all,tabJoin,by='fish_id')

metrics.fast_slow.plot <- left_join(metrics.fast_slow.all,tabJoin,by='fish_id')

windows()
ggplot(subset(metrics.model.plot,ramping == 'fast' & freqChan == 'all'),aes(x=as.factor(angleMatch),y=cor))+
  geom_boxplot()+
  ylim(0,10)

windows()
ggplot(subset(metrics.model.plot,ramping == 'slow' & freqChan == '70khz' & angleMatch == 80),aes(x=fish_length,y=mean))+
  geom_point()

windows()
ggplot(subset(metrics.fast_slow.plot,freqChan == 'all' & fish_id == 'P03'),aes(x=angleMatch,y=cor))+
  geom_line(data=subset(metrics.fast_slow.plot,freqChan == 'all' & fish_id == 'P03'),aes(x=angleMatch,y=cor),col='red')

windows()
ggplot()+
  geom_line(data=subset(metrics.fast_slow.plot,freqChan == '70khz' & fish_id == 'S08'),aes(x=angleMatch,y=mean),col='blue')+
  geom_line(data=subset(metrics.fast_slow.plot,freqChan == '70khz' & fish_id == 'P01'),aes(x=angleMatch,y=mean),col='red')+
  geom_line(data=subset(metrics.fast_slow.plot,freqChan == '70khz' & fish_id == 'S27'),aes(x=angleMatch,y=mean),col='yellow')

################################################################################
## Dump - old
################################################################################

for(idx in 1:dim(summaryTab)[1]){
  KRMCurrent <- read.csv(file.path(myPath$resultsPath,'KRM',
                                   paste0('krm_',summaryTab$fish_id[idx],'.csv')))

  measuName <- paste0(summaryTab$fish_id[idx],'_',
                      summaryTab$mode[idx],'_',
                      summaryTab$ramping[idx],'_',
                      summaryTab$orientation[idx],'_',
                      summaryTab$resolution[idx])
  
  #print(measuName)
  
  ###################################
  # load 200khz results
  ###################################
  measuName.200khz  <- paste0(measuName,'_200khz')
  TS.cor.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                      paste0(measuName.200khz,'_TSCorr.csv')),check.names=FALSE)
  TS.cor.200khz$angleMatch <- TS.cor.200khz$angle-270
  TS.cor.200khz$angleMatch[TS.cor.200khz$angleMatch < 0] <- TS.cor.200khz$angleMatch[TS.cor.200khz$angleMatch < 0]+360
  TS.cor.200khz$freqChan <- '200khz'
  
  TS.kde.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                      paste0(measuName.200khz,'_TSkde.csv')),check.names=FALSE)
  TS.kde.200khz <- TS.kde.200khz %>% pivot_longer(!angle & !TSmesh,names_to = "frequency", values_to = "density")
  TS.kde.200khz$frequency <- as.numeric(TS.kde.200khz$frequency)
  TS.kde.200khz$angleMatch <- TS.kde.200khz$angle-270
  TS.kde.200khz$angleMatch[TS.kde.200khz$angleMatch < 0] <- TS.kde.200khz$angleMatch[TS.kde.200khz$angleMatch < 0]+360
  TS.kde.200khz$freqChan <- '200khz'
  
  TS.stats.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                        paste0(measuName.200khz,'_TSStats.csv')),check.names=FALSE)
  TS.stats.200khz <- TS.stats.200khz %>% pivot_longer(!angle & !percentiles,names_to = "frequency", values_to = "TS")
  TS.stats.200khz$frequency <- as.numeric(TS.stats.200khz$frequency)
  TS.stats.200khz$angleMatch <- TS.stats.200khz$angle-270
  TS.stats.200khz$angleMatch[TS.stats.200khz$angleMatch < 0] <- TS.stats.200khz$angleMatch[TS.stats.200khz$angleMatch < 0]+360
  TS.stats.200khz$freqChan <- '200khz'
  
  ###################################
  # load 70khz results
  ###################################
  measuName.70khz   <- paste0(measuName,'_70khz')
  TS.cor.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                     paste0(measuName.70khz,'_TSCorr.csv')),check.names=FALSE)
  TS.cor.70khz$angleMatch <- TS.cor.70khz$angle-270
  TS.cor.70khz$angleMatch[TS.cor.70khz$angleMatch < 0] <- TS.cor.70khz$angleMatch[TS.cor.70khz$angleMatch < 0]+360
  TS.cor.70khz$freqChan <- '70khz'
  
  TS.kde.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                     paste0(measuName.70khz,'_TSkde.csv')),check.names=FALSE)
  TS.kde.70khz <- TS.kde.70khz %>% pivot_longer(!angle & !TSmesh,names_to = "frequency", values_to = "density")
  TS.kde.70khz$frequency <- as.numeric(TS.kde.70khz$frequency)
  TS.kde.70khz$angleMatch <- TS.kde.70khz$angle-270
  TS.kde.70khz$angleMatch[TS.kde.70khz$angleMatch < 0] <- TS.kde.70khz$angleMatch[TS.kde.70khz$angleMatch < 0]+360
  TS.kde.70khz$freqChan <- '70khz'
  
  TS.stats.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                       paste0(measuName.70khz,'_TSStats.csv')),check.names=FALSE)
  TS.stats.70khz <- TS.stats.70khz %>% pivot_longer(!angle & !percentiles,names_to = "frequency", values_to = "TS")
  TS.stats.70khz$frequency <- as.numeric(TS.stats.70khz$frequency)
  TS.stats.70khz$angleMatch <- TS.stats.70khz$angle-270
  TS.stats.70khz$angleMatch[TS.stats.70khz$angleMatch < 0] <- TS.stats.70khz$angleMatch[TS.stats.70khz$angleMatch < 0]+360
  TS.stats.70khz$freqChan <- '70khz'
  
  ###################################
  # combine tables
  ###################################
  TS.cor.all <- rbind(TS.cor.70khz,TS.cor.200khz)
  TS.stats.all <- rbind(TS.stats.70khz,TS.stats.200khz)
  TS.kde.all <- rbind(TS.kde.70khz,TS.kde.200khz)
  
  uniqueAngles <- unique(TS.cor.all$angleMatch)
  
  ###################################
  # plotting
  ###################################
  dir.create(file.path(myPath$figuresPath,measuName),showWarnings = FALSE)
  
  pdf(file.path(myPath$figuresPath,measuName,paste0(measuName,'TS correlation.pdf')))
  
  ggplot(TS.cor.all,aes(x=angleMatch,y=TSCorr,col=freqChan))+
    geom_line()+
    ggtitle('TS(f) correlation within angle bin')+
    theme_bw()
  
  dev.off()
  
  # plotting against TS density
  pdf(file.path(myPath$figuresPath,measuName,paste0(measuName,'TS kde.pdf')))
  for(angleSelect in uniqueAngles){
    plot.200khz <- subset(TS.kde.200khz,angleMatch == angleSelect)
    plot.200khz$density[plot.200khz$density < 1e-2] <- NA
    plot.200khz$density <-  plot.200khz$density/max(plot.200khz$density,na.rm=TRUE)
    plot.70khz  <- subset(TS.kde.70khz,angleMatch == angleSelect)
    plot.70khz$density[plot.70khz$density < 1e-2] <- NA
    plot.70khz$density <-  plot.70khz$density/max(plot.70khz$density,na.rm=TRUE)
    
    if(angleSelect >= 65 & angleSelect <= 115){
      plot.krm <- subset(KRMCurrent,theta == angleSelect)
      
      windows()
      ggplot()+
        geom_line(data=plot.krm,
                  aes(x=frequency/1e3,y=TS),size=1,col='black')
      
      print(ggplot()+
        geom_tile(data=plot.200khz,
                  aes(x=frequency/1e3,y=TSmesh,fill=density))+
        geom_tile(data=plot.70khz,
                  aes(x=frequency/1e3,y=TSmesh,fill=density))+
        geom_line(data=plot.krm,
                  aes(x=frequency/1e3,y=TSmod),size=1,col='black')+
        scale_fill_viridis_b()+
        xlim(50,250)+
        ylim(-70,-20)+
        ggtitle(paste0(measuName,'<<-->> angle:',angleSelect))+
        theme_bw())
    }else{
      print(ggplot()+
        geom_tile(data=plot.200khz,
                  aes(x=frequency/1e3,y=TSmesh,fill=density))+
        geom_tile(data=plot.70khz,
                  aes(x=frequency/1e3,y=TSmesh,fill=density))+
        scale_fill_viridis_b()+
        xlim(50,250)+
        ylim(-70,-20)+
        ggtitle(paste0(measuName,'<<-->> angle:',angleSelect))+
        theme_bw())
    }
  }
  dev.off()
  
  # plotting against TS density
  pdf(file.path(myPath$figuresPath,measuName,paste0(measuName,'TS stats.pdf')))
  for(angleSelect in uniqueAngles){
    plot.all <- subset(TS.stats.all,angleMatch == angleSelect)
    plot.all$percentiles <- as.character(plot.all$percentiles)
    plot.all <- plot.all %>% pivot_wider(values_from=TS,names_from=percentiles)
    
    if(angleSelect >= 65 & angleSelect <= 115){
      plot.krm <- subset(KRMCurrent,theta == angleSelect)
      
      print(ggplot()+
        geom_ribbon(data=plot.all,
                    aes(x=frequency/1e3,ymin=`25`,ymax=`75`,fill=freqChan))+
        geom_line(data=plot.krm,
                  aes(x=frequency/1e3,y=TS),size=1,col='black')+
        ggtitle(paste0(measuName,'<<-->> angle:',angleSelect))+
        xlim(50,250)+
        ylim(-70,-20)+
        theme_bw())
    }else{
      print(ggplot()+
        geom_ribbon(data=plot.all,
                    aes(x=frequency/1e3,ymin=`25`,ymax=`75`,fill=freqChan))+
        ggtitle(paste0(measuName,'<<-->> angle:',angleSelect))+
        xlim(50,250)+
        ylim(-70,-20)+
        theme_bw())
    }
  }
  dev.off()
}