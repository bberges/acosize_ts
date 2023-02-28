plotFun <- function(myPath,summaryTab,idx){
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
      
      print(ggplot()+
              geom_tile(data=plot.200khz,
                        aes(x=frequency/1e3,y=TSmesh,fill=density))+
              geom_tile(data=plot.70khz,
                        aes(x=frequency/1e3,y=TSmesh,fill=density))+
              geom_line(data=plot.krm,
                        aes(x=frequency/1e3,y=TS),size=1,col='black')+
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