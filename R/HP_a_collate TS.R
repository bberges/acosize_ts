library(ggplot2)
library(tidyverse)
library(ggbiplot)
library(mgcv)
library(ale)

rm(list = ls())

path <- 'G:/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     rPath=file.path(".","R"))

freq <- 70
processRaw <- T

#--------------------------
# read/create data
#--------------------------
if(processRaw){
  files.spectra <- list.files(path = file.path(myPath$resultsPath,
                                       'measurements',
                                       paste0(freq,'khz')),
                      pattern = 'TSMat.csv')

  files.time <- list.files(path = file.path(myPath$resultsPath,
                                               'measurements',
                                               paste0(freq,'khz')),
                              pattern = 'TimeMat.csv')

  flagFirst <- T
  for(file in files.spectra){
    TS <- read.csv(file.path(myPath$resultsPath,
                             'measurements',
                             paste0(freq,'khz'),
                             file),check.names = F)

    if(flagFirst){
      TS.all <- TS
      flagFirst <- F
    }else{
      TS.all <- rbind(TS.all,TS)
    }
  }

  flagFirst <- T
  for(file in files.time){
    TF <- read.csv(file.path(myPath$resultsPath,
                             'measurements',
                             paste0(freq,'khz'),
                             file),check.names = F)

    if(flagFirst){
      TF.all <- TF
      flagFirst <- F
    }else{
      TF.all <- rbind(TF.all,TF)
    }
  }

  data.all <- left_join(TS.all,TF.all,by=colnames(TS.all)[colnames(TS.all) %in% colnames(TF.all)])

  data.all$angleGroup <- cut(data.all$angle,
                           seq(min(data.all$angle),max(data.all$angle)+1,by=3),
                           include.lowest = T)
  data.all <- data.all %>% relocate(angleGroup)
  data.all$ID <- paste0(data.all$fish_id,'-',
                                data.all$ramping,'-',
                                data.all$orientation,'-',
                                data.all$angle,'-',
                                data.all$spectra_id)
  rownames(data.all) <- data.all$ID

  freq.col <- !is.na(as.numeric(colnames(data.all)))


  for(idxRamping in unique(data.all$ramping)){
    print(idxRamping)
    for(idxOrient in unique(data.all$orientation)){
      print(idxOrient)
      # filtering outliers
      for(idxID in unique(data.all$fish_id)){
        print(idxID)
        for(idxGroup in unique(data.all$angleGroup)){
          #print(idxGroup)
          idx.sel <- data.all$angleGroup == idxGroup &
            data.all$fish_id == idxID &
            data.all$ramping == idxRamping &
            data.all$orientation == idxOrient

          #TS.sub <- TS.all[idx.sel,]
          dist.mat <- (dist(data.all[idx.sel,freq.col], method = "euclidean"))
          df <- reshape2::melt(as.matrix(dist.mat), varnames = c("row", "col"))
          df <- df %>% group_by(row) %>% summarize(value=mean(value))

          # TS.sub <- TS.all[idx.sel,]
          # TS.sub <- TS.sub %>% pivot_longer(cols = colnames(TS.sub[,13:dim(TS.sub)[2]]),values_to = 'TS',names_to = 'freq')
          # TS.sub$freq <- as.numeric(TS.sub$freq)
          #
          # windows()
          # ggplot(TS.sub,aes(x=freq,y=TS,col=as.factor(spectra_id)))+
          #   geom_line()

          upper_bound <- median(df$value) + 3 * mad(df$value, constant = 1)
          idx.filt <- which(df$value > upper_bound)

          # TS.sub <- TS.sub[!(TS.sub$spectra_id %in% df$row[idx.filt]),]
          #
          # windows()
          # ggplot(TS.sub,aes(x=freq,y=TS,col=as.factor(spectra_id)))+
          #   geom_line()

          data.all <- data.all[!(data.all$ID %in% df$row[idx.filt]),]
        }
      }
    }
  }

  save(data.all,file = file.path(myPath$resultsPath,paste0('data_',freq,'khz.RData')))
}else{
  load(file = file.path(myPath$resultsPath,paste0('data_',freq,'khz.RData')))
}
