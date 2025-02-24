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

source(file.path(myPath$rPath,'utilities_calc.R'))

summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]
summaryTab <- summaryTab[summaryTab$orientation == 'ventral',]

uniqueAnimals <- unique(summaryTab$fish_id)

distanceMeasure <- 'eucl'
windowing <- 'variableWindow' # variableWindow fixedWindow
dir.create(file.path(myPath$figuresPath,paste0(windowing,'_',distanceMeasure)),showWarnings = FALSE)
myPath$figuresPath <- file.path(myPath$figuresPath,paste0(windowing,'_',distanceMeasure))

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
                        unique(currentTab$resolution),'_',
                        windowing)

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

    TS.mat.200khz <- read.csv(file.path(myPath$resultsPath,'measurements','200khz',
                                          paste0(measuName.200khz,'_TSMat.csv')),check.names=FALSE)
    TS.mat.200khz$spectra_id.str.angle <- paste0(TS.mat.200khz$angle,'_',idxRamping,'_',TS.mat.200khz$spectra_id)
    TS.mat.200khz <- TS.mat.200khz %>% pivot_longer(!angle & !spectra_id & !spectra_id.str.angle,names_to = "frequency", values_to = "TS")
    TS.mat.200khz$frequency <- as.numeric(TS.mat.200khz$frequency)
    TS.mat.200khz$angleMatch <- TS.mat.200khz$angle-270
    TS.mat.200khz$angleMatch[TS.mat.200khz$angleMatch < 0] <- TS.mat.200khz$angleMatch[TS.mat.200khz$angleMatch < 0]+360
    TS.mat.200khz$spectra_id.str <- paste0(TS.mat.200khz$angleMatch,'_',idxRamping,'_',TS.mat.200khz$spectra_id)
    TS.mat.200khz$freqChan <- '200khz'
    TS.mat.200khz$ramping <- idxRamping
    TS.mat.200khz <- TS.mat.200khz %>% select(-c('spectra_id.str.angle'))

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

    TS.mat.70khz <- read.csv(file.path(myPath$resultsPath,'measurements','70khz',
                                        paste0(measuName.70khz,'_TSMat.csv')),check.names=FALSE)
    TS.mat.70khz$spectra_id.str.angle <- paste0(TS.mat.70khz$angle,'_',idxRamping,'_',TS.mat.70khz$spectra_id)
    TS.mat.70khz <- TS.mat.70khz %>% pivot_longer(!angle & !spectra_id & !spectra_id.str.angle,names_to = "frequency", values_to = "TS")
    TS.mat.70khz$frequency <- as.numeric(TS.mat.70khz$frequency)
    TS.mat.70khz$angleMatch <- TS.mat.70khz$angle-270
    TS.mat.70khz$angleMatch[TS.mat.70khz$angleMatch < 0] <- TS.mat.70khz$angleMatch[TS.mat.70khz$angleMatch < 0]+360
    TS.mat.70khz$spectra_id.str <- paste0(TS.mat.70khz$angleMatch,'_',idxRamping,'_',TS.mat.70khz$spectra_id)
    TS.mat.70khz$freqChan <- '70khz'
    TS.mat.70khz$ramping <- idxRamping
    TS.mat.70khz <- TS.mat.70khz %>% select(-c('spectra_id.str.angle'))

    ###################################
    # combine tables
    ###################################
    if(flagFirst){
      TS.cor.all <- rbind(TS.cor.70khz,TS.cor.200khz)
      TS.stats.all <- rbind(TS.stats.70khz,TS.stats.200khz)
      TS.kde.all <- rbind(TS.kde.70khz,TS.kde.200khz)
      TS.mat.all <- rbind(TS.mat.70khz,TS.mat.200khz)
      flagFirst <- FALSE
    }else{
      TS.cor.all <- rbind(TS.cor.all,TS.cor.70khz,TS.cor.200khz)
      TS.stats.all <- rbind(TS.stats.all,TS.stats.70khz,TS.stats.200khz)
      TS.kde.all <- rbind(TS.kde.all,TS.kde.70khz,TS.kde.200khz)
      TS.mat.all <- rbind(TS.mat.all,TS.mat.70khz,TS.mat.200khz)
    }
  }

  # temp <- subset(TS.mat.all,angleMatch==0 & ramping=='slow')
  # temp$TS[temp$frequency == 87000] <- NA
  #
  # scaling_factor <- 1
  # png(file.path(myPath$figuresPath,'spectra_ind_0.png'),
  #     width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)
  #
  # p <- ggplot(temp,aes(x=frequency*1e-3,y=TS,col=as.factor(spectra_id)))+
  #       geom_line()+
  #       theme_bw()+
  #       theme(legend.position = 'none')+
  #       xlab('Frequency (kHz)')+
  #       ylab(expression('TS (dB re 1m'^'2'*')'))+
  #       ylim(-70,-30)
  #
  # print(p)
  # dev.off()
  #
  # out$metrics.slow_fast$idxBestFast

  #angleSel <- 115
  #windows()
  #ggplot()+
  #  geom_tile(data=subset(TS.kde.all,ramping=='fast' & angleMatch == angleSel & freqChan =='70khz'),
  #            aes(x=frequency,y=TSmesh,fill=density))+
  #  geom_tile(data=subset(TS.kde.all,ramping=='fast' & angleMatch == angleSel & freqChan =='200khz'),
  #            aes(x=frequency,y=TSmesh,fill=density))+
  #  geom_line(data=subset(TS.mat.all,ramping=='fast' & angleMatch == angleSel & freqChan =='70khz'),aes(x=frequency,y=TS,col=as.factor(spectra_id)),alpha=0.3)+
  #  geom_line(data=subset(TS.mat.all,ramping=='fast' & angleMatch == angleSel & freqChan =='200khz'),aes(x=frequency,y=TS,col=as.factor(spectra_id)),alpha=0.3)+
  #  geom_line(data=subset(KRMCurrent,theta==angleSel),aes(x=frequency,y=TS),linewidth=1)+
  #  scale_fill_viridis_b()+
  #  ylim(-70,-25)+
  #  xlim(50*1e3,280*1e3)

  KRMCurrent.filt <- KRMCurrent[KRMCurrent$frequency %in% unique(TS.stats.all$frequency),]
  uniqueAngles.model <- unique(KRMCurrent.filt$theta)

  freq.70khz <- unique(TS.stats.all$frequency[TS.stats.all$freqChan == '70khz'])
  freq.200khz <- unique(TS.stats.all$frequency[TS.stats.all$freqChan == '200khz'])

  idx70khz <- unique(TS.stats.all$frequency) %in% freq.70khz
  idx200khz <- unique(TS.stats.all$frequency) %in% freq.200khz

  TS.mat.wide <- TS.mat.all  %>% group_by(angle,frequency,angleMatch,ramping,spectra_id.str) %>% summarize(TS=TS) %>% pivot_wider(names_from=frequency,values_from=TS)
  uniqueAngles <- unique(TS.mat.wide$angleMatch)

  # open pdf for plotting
  pdf(file.path(myPath$figuresPath,paste0(measuName,'_TS spectra.pdf')))

  flagFirst <- TRUE
  for(idxAngle in uniqueAngles){
    print(idxAngle)
    TS.mat.wide.current <- subset(TS.mat.wide,angleMatch == idxAngle)
    TS.mat.wide.current.info <- TS.mat.wide.current %>% ungroup() %>% select(c('ramping','angle','angleMatch','spectra_id.str'))
    TS.mat.wide.current <- TS.mat.wide.current %>% ungroup() %>% select(-c('angle','angleMatch','spectra_id.str'))

    KRMCurrent.current <- subset(KRMCurrent.filt,theta == idxAngle) %>% dplyr::select(c('frequency','TS'))# %>% pivot_longer(names_from=frequency,values_from=TS)

    angleMatch <- idxAngle

    if(distanceMeasure == 'dtw'){
      out <- computeDistDTWMat(idxAngle,
                               KRMCurrent.current,
                               TS.mat.wide.current,
                               TS.mat.wide.current.info,
                               idx70khz,idx200khz)
    }else if(distanceMeasure == 'eucl'){
      out <- computeDistEuclMat(idxAngle,
                                KRMCurrent.current,
                                TS.mat.wide.current,
                                TS.mat.wide.current.info,
                                idx70khz,idx200khz)
    }

    # # pull out matching spectra
    flagFirstSpectra <- TRUE
    for(idxFreq in c('70khz','200khz')){ # 'all'
      if(idxFreq == '70khz'){
        freq.filt <- freq.70khz
        freq.sel <- '70khz'
      }else if(idxFreq == '200khz'){
        freq.filt <- freq.200khz
        freq.sel <- '200khz'
      }else if(idxFreq =='all'){
        freq.filt <- c(freq.70khz,freq.200khz)
        freq.sel <- c('70khz','200khz')
      }

      if(dim(subset(out$metrics.slow_fast,freqChan == idxFreq))[1] != 0){
        TS.best.fast <- subset(TS.mat.all,freqChan %in% freq.sel &
                                 spectra_id.str == subset(out$metrics.slow_fast,freqChan == idxFreq)$idxBestFast)
        TS.best.slow <- subset(TS.mat.all,freqChan %in% freq.sel &
                                 spectra_id.str == subset(out$metrics.slow_fast,freqChan == idxFreq)$idxBestSlow)
        spectra.best <- rbind(TS.best.fast,TS.best.slow)
        spectra.best <- subset(spectra.best,frequency %in% freq.filt)
        spectra.best$cat <- 'slow_fast'

        if(flagFirstSpectra){
          spectra.best.all <- spectra.best
          spectra.best.all$freqChan <- factor(spectra.best.all$freqChan,levels=c('70khz','200khz'))
          spectra.best.all$cat <- factor(spectra.best.all$cat,levels=c('model','slow_fast'))
          spectra.best.all$ramping <- factor(spectra.best.all$ramping,levels=c('slow','fast','model'))
          flagFirstSpectra <- FALSE
        }else{
          spectra.best.all <- rbind(spectra.best.all,spectra.best)
        }
      }

      if(dim(subset(out$metrics.model,freqChan == idxFreq))[1] != 0){
        metrics.model.current <- subset(out$metrics.model,freqChan == idxFreq)

        for(idxSpecRamp in unique(metrics.model.current$ramping)){
          str.spectra.best <- subset(metrics.model.current,ramping == idxSpecRamp)$idxBest
          spectra.best <- subset(TS.mat.all,freqChan %in% freq.sel &
                                   spectra_id.str == str.spectra.best)
          spectra.best <- subset(spectra.best,frequency %in% freq.filt)
          spectra.best$cat <- 'model'

          if(flagFirstSpectra){
            spectra.best.all <- spectra.best
            spectra.best.all$freqChan <- factor(spectra.best.all$freqChan,levels=c('70khz','200khz'))
            spectra.best.all$cat <- factor(spectra.best.all$cat,levels=c('model','slow_fast'))
            spectra.best.all$ramping <- factor(spectra.best.all$ramping,levels=c('slow','fast','model'))
            flagFirstSpectra <- FALSE
          }else{
            spectra.best.all <- rbind(spectra.best.all,spectra.best)
          }
        }

        spectrum.model <- KRMCurrent.current
        spectrum.model$angle <- unique(subset(TS.mat.wide,angleMatch == idxAngle)$angle)
        spectrum.model$spectra_id <- -1
        spectrum.model$angleMatch <- idxAngle
        spectrum.model$freqChan   <- NA
        spectrum.model$freqChan[spectrum.model$frequency %in% freq.70khz] <- '70khz'
        spectrum.model$freqChan[spectrum.model$frequency %in% freq.200khz] <- '200khz'
        spectrum.model$spectra_id.str <- paste0(idxAngle,'_model_0')
        spectrum.model$ramping <- 'model'
        spectrum.model$cat <- 'model'
        spectrum.model <- subset(spectrum.model,freqChan == idxFreq)

        spectra.best.all <- rbind(spectra.best.all,spectrum.model)
      }
    }

    if(exists('spectra.best.all')){
      print(ggplot(spectra.best.all,aes(x=frequency*1e-3,y=TS,col=ramping))+
              geom_line()+
              facet_grid(cat~freqChan,scales = 'free_x',drop=FALSE)+
              labs(title=paste0(idxAnimal,' - angle=',idxAngle),
                   x='Frequency (kHz)',
                   y='TS (dB)'))
              #xlim(50,280)+
              #ylim(-75,-25))
      rm(spectra.best.all)
    }


    if(flagFirst){
      metrics.slow_fast <- out$metrics.slow_fast
      metrics.model     <- out$metrics.model
      flagFirst <- FALSE
    }else{
      metrics.slow_fast <- rbind(metrics.slow_fast,out$metrics.slow_fast)
      metrics.model     <- rbind(metrics.model,out$metrics.model)
    }
  }

  # closing pdf for plotting
  dev.off()


  # ggplot(subset(metrics.slow_fast,freqChan =='70khz'),aes(x=angleMatch,y=log(meanDist)))+
  #   geom_point()
  #
  # ggplot()+
  #   geom_point(data=subset(metrics.model,freqChan =='200khz' & ramping == 'slow'),aes(x=angleMatch,y=log(sdDist)),col='red')+
  #   geom_point(data=subset(metrics.model,freqChan =='70khz' & ramping == 'slow'),aes(x=angleMatch,y=log(sdDist)),col='blue')+
  #   xlim(65,115)
  #
  # ggplot(subset(metrics.model,freqChan =='70khz' & ramping == 'fast'),aes(x=angleMatch,y=(meanDistTrunc25)))+
  #   geom_line()
  #
  # ggplot(subset(metrics.slow_fast,freqChan =='200khz'),aes(x=angleMatch,y=(sdDist)))+
  #   geom_point()
  #
  #
  # ggplot(subset(metrics.model.dB,freqChan =='70khz' & ramping == 'slow'),aes(x=angleMatch,y=log(minDist)))+
  #   geom_point()

  ##################################################################
  # plotting individual spectra best fast/slow match
  ##################################################################
  metrics.slow_fast.alt <- metrics.slow_fast[!is.na(metrics.slow_fast$meanDist),]
  uniqueAngles <- unique(metrics.slow_fast$angleMatch)

  idxFreq <- '200khz'
  flagFirstTemp <- TRUE
  for(idxAngle in uniqueAngles){
    metrics.slow_fast.filt <- subset(metrics.slow_fast.alt,angleMatch == idxAngle & freqChan == idxFreq)

    if(dim(metrics.slow_fast.filt)[1] != 0){
      mySpectra.slow <- subset(TS.mat.all,freqChan == idxFreq & spectra_id.str == metrics.slow_fast.filt$idxBestSlow)
      mySpectra.fast <- subset(TS.mat.all,freqChan == idxFreq & spectra_id.str == metrics.slow_fast.filt$idxBestFast)

      mySpectra <- rbind(mySpectra.slow,mySpectra.fast)

      if(flagFirstTemp){
        spectra.plot <- mySpectra
        flagFirstTemp <- FALSE
      }else{
        spectra.plot <- rbind(spectra.plot,mySpectra)
      }
    }
  }

  scaling_factor <- 2
  png(file.path(myPath$figuresPath,paste0('TS_measurements_',idxAnimal,'_200khz.png')),
      width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

  p <- ggplot(data=spectra.plot, aes(x=frequency/1e3, y=angleMatch, fill=TS))+
        geom_tile()+
        scale_fill_gradientn(colors=rev(pals::brewer.spectral(15)),
                             limits=c(-60,-30), oob=scales::squish,
                             name=expression(TS~"(dB re 1"*m^2*")"))+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        ylab(expression(theta~"(°)"))+
        xlab("Frequency (kHz)")+
        theme_classic()+
        theme(legend.position = "bottom",
              text=element_text(size=16),
              legend.key.width = unit(2, 'cm'))+
        facet_wrap(~ramping)

  print(p)
  dev.off()


  idxFreq <- '70khz'
  flagFirstTemp <- TRUE
  for(idxAngle in uniqueAngles){
    metrics.slow_fast.filt <- subset(metrics.slow_fast.alt,angleMatch == idxAngle & freqChan == idxFreq)

    if(dim(metrics.slow_fast.filt)[1] != 0){
      mySpectra.slow <- subset(TS.mat.all,freqChan == idxFreq & spectra_id.str == metrics.slow_fast.filt$idxBestSlow)
      mySpectra.fast <- subset(TS.mat.all,freqChan == idxFreq & spectra_id.str == metrics.slow_fast.filt$idxBestFast)

      mySpectra <- rbind(mySpectra.slow,mySpectra.fast)

      if(flagFirstTemp){
        spectra.plot <- mySpectra
        flagFirstTemp <- FALSE
      }else{
        spectra.plot <- rbind(spectra.plot,mySpectra)
      }
    }
  }

  scaling_factor <- 2
  png(file.path(myPath$figuresPath,paste0('TS_measurements_',idxAnimal,'_70khz.png')),
      width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

  p <- ggplot(data=spectra.plot, aes(x=frequency/1e3, y=angleMatch, fill=TS))+
    geom_tile()+
    scale_fill_gradientn(colors=rev(pals::brewer.spectral(15)),
                         limits=c(-60,-30), oob=scales::squish,
                         name=expression(TS~"(dB re 1"*m^2*")"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    ylab(expression(theta~"(°)"))+
    xlab("Frequency (kHz)")+
    theme_classic()+
    theme(legend.position = "bottom",
          text=element_text(size=16),
          legend.key.width = unit(2, 'cm'))+
    facet_wrap(~ramping)

  print(p)
  dev.off()

  ##################################################################
  # plotting individual spectra against model
  ##################################################################
  metrics.model.alt <- metrics.model[!is.na(metrics.model$meanTS.model),]
  uniqueAngles <- unique(metrics.model$angleMatch)


  scaling_factor <- 2
  png(file.path(myPath$figuresPath,paste0('TS_model vs measurements_',idxAnimal,'_70khz.png')),
      width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

  idxFreq <- '70khz'
  idxRamping  <- 'fast'
  flagFirstTemp <- TRUE
  for(idxAngle in uniqueAngles){
    metrics.model.filt <- subset(metrics.model.alt,angleMatch == idxAngle & freqChan == idxFreq & ramping == idxRamping)

    if(dim(metrics.model.filt)[1] != 0){
      mySpectra <- subset(TS.mat.all,freqChan == idxFreq & ramping == idxRamping & spectra_id.str == metrics.model.filt$idxBest)

      if(flagFirstTemp){
        spectra.plot <- mySpectra
        flagFirstTemp <- FALSE
      }else{
        spectra.plot <- rbind(spectra.plot,mySpectra)
      }
    }
  }

  spectra.plot <- spectra.plot %>% select(c('frequency','angleMatch','TS'))
  spectra.plot$type <- 'measurements'
  KRMCurrent.plot <- KRMCurrent.filt[KRMCurrent.filt$frequency %in% freq.70khz,]
  KRMCurrent.plot$angleMatch <- KRMCurrent.plot$theta
  KRMCurrent.plot <- KRMCurrent.plot %>% select(c('frequency','angleMatch','TS'))
  KRMCurrent.plot$type <- 'model'
  spectra.plot <- rbind(spectra.plot,KRMCurrent.plot)

  p <- ggplot(data=spectra.plot, aes(x=frequency/1e3, y=angleMatch, fill=TS))+
        geom_tile()+
        scale_fill_gradientn(colors=rev(pals::brewer.spectral(15)),
                             limits=c(-60,-30), oob=scales::squish,
                             name=expression(TS~"(dB re 1"*m^2*")"))+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        ylab(expression(theta~"(°)"))+
        xlab("Frequency (kHz)")+
        theme_classic()+
        theme(legend.position = "bottom",
              text=element_text(size=16),
              legend.key.width = unit(2, 'cm'))+
    facet_wrap(~type)

  print(p)
  dev.off()


  scaling_factor <- 2
  png(file.path(myPath$figuresPath,paste0('TS_model vs measurements_',idxAnimal,'_200khz.png')),
      width = 12*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

  idxFreq <- '200khz'
  idxRamping  <- 'fast'
  flagFirstTemp <- TRUE
  for(idxAngle in uniqueAngles){
    metrics.model.filt <- subset(metrics.model.alt,angleMatch == idxAngle & freqChan == idxFreq & ramping == idxRamping)

    if(dim(metrics.model.filt)[1] != 0){
      mySpectra <- subset(TS.mat.all,freqChan == idxFreq & ramping == idxRamping & spectra_id.str == metrics.model.filt$idxBest)

      if(flagFirstTemp){
        spectra.plot <- mySpectra
        flagFirstTemp <- FALSE
      }else{
        spectra.plot <- rbind(spectra.plot,mySpectra)
      }
    }
  }

  spectra.plot <- spectra.plot %>% select(c('frequency','angleMatch','TS'))
  spectra.plot$type <- 'measurements'
  KRMCurrent.plot <- KRMCurrent.filt[KRMCurrent.filt$frequency %in% freq.200khz,]
  KRMCurrent.plot$angleMatch <- KRMCurrent.plot$theta
  KRMCurrent.plot <- KRMCurrent.plot %>% select(c('frequency','angleMatch','TS'))
  KRMCurrent.plot$type <- 'model'
  spectra.plot <- rbind(spectra.plot,KRMCurrent.plot)

  p <- ggplot(data=spectra.plot, aes(x=frequency/1e3, y=angleMatch, fill=TS))+
    geom_tile()+
    scale_fill_gradientn(colors=rev(pals::brewer.spectral(15)),
                         limits=c(-60,-30), oob=scales::squish,
                         name=expression(TS~"(dB re 1"*m^2*")"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    ylab(expression(theta~"(°)"))+
    xlab("Frequency (kHz)")+
    theme_classic()+
    theme(legend.position = "bottom",
          text=element_text(size=16),
          legend.key.width = unit(2, 'cm'))+
  facet_wrap(~type)

  print(p)
  dev.off()


  #####################################################
  metrics.slow_fast$fish_id <- idxAnimal
  metrics.model$fish_id     <- idxAnimal

  # windows()
  # ggplot()+
  #  geom_line(data=metrics.slow_fast,aes(x=angleMatch,y=meanDistTrunc25,col=freqChan))
  #
  # windows()
  # ggplot()+
  #  geom_line(data=subset(metrics.model,!is.na(minDist)),aes(x=angleMatch,y=meanDist,col=freqChan))+
  #  facet_wrap(~ramping)

  if(flagFirstAnimal){
    metrics.slow_fast.all <- metrics.slow_fast
    metrics.model.all     <- metrics.model
    flagFirstAnimal       <- FALSE
  }else{
    metrics.slow_fast.all <- rbind(metrics.slow_fast.all,metrics.slow_fast)
    metrics.model.all     <- rbind(metrics.model.all,metrics.model)
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

save(metrics.slow_fast.all,
     metrics.model.all,
     file = file.path(myPath$resultsPath,paste0('acosize_TS_summary results_',
                                                windowing, '_',
                                                distanceMeasure,'.RData')))
