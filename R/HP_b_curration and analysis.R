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

freq <- 200
orientation.sel <- 'ventral'
ramping.sel <- 'slow'
type.sel <- c('SAI')

load(file = file.path(myPath$resultsPath,paste0('data_',freq,'khz.RData')))

freq.col <- !is.na(as.numeric(colnames(data.all)))

#--------------------------
# build data set
#--------------------------
data.sub <- subset(data.all,ramping == ramping.sel &
                     orientation == orientation.sel &
                     type %in% type.sel)

data.mean <- apply(data.sub[,freq.col],1,'mean')
data.norm <- data.sub
data.norm[,freq.col] <- sweep(sweep(data.norm[,freq.col],
                                    1,apply(data.norm[,freq.col],1,'mean'),
                                    '-'),
                              1,apply(data.norm[,freq.col],1,'sd'),
                              '/')

data.PCA <- prcomp(data.norm[,freq.col], center = TRUE,scale. = TRUE)

df <- data.frame(angle = data.sub$angle,
                 length = data.sub$fish_length,
                 ID = as.factor(data.sub$fish_id),
                 SNR = data.sub$SNR_sig,
                 meanVal = apply(data.sub[,freq.col],1,'mean'),
                 sdVal = apply(data.sub[,freq.col],1,'sd'),
                 PC1 = data.PCA$x[,1],
                 PC2 = data.PCA$x[,2],
                 nullsN = data.sub$nullsN,
                 nullsProm = data.sub$nullsProm,
                 nullsWidthMean = data.sub$nullsWidthMean,
                 mu = data.sub$mu,
                 sigma = data.sub$sigma,
                 skew = data.sub$skew,
                 kurt = data.sub$kurt,
                 Npeaks=data.sub$Npeaks,
                 widthPeaks=data.sub$widthPeaks,
                 promPeaks=data.sub$promPeaks)

df$lengthGroup <- cut(df$length,
                      seq(25,55,by=5),
                      include.lowest = T)

df$lengthGroupNum <- as.numeric(df$lengthGroup)

df$angleGroup <- cut(df$angle,
                     seq(min(df$angle),360,by=36),
                     include.lowest = T)

df <- df %>% filter(!if_any(everything(), is.na))

write.csv(file = file.path(myPath$resultsPath,
                           paste0('acosize_',
                                  freq,
                                  '_',
                                  type.sel,
                                  '_features and targets.csv')),
          df)


#----------------------------------------------------------------------
#----------------------------------------------------------------------


windows()
ggplot(df,aes(x=angle,y=nullsN,group=angleGroup))+
  stat_summary(fun.y = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "pointrange")

windows()
ggplot(df,aes(x=length,y=nullsN,group=lengthGroup))+
  stat_summary(fun.y = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "pointrange")+
  facet_wrap(~angleGroup,scales='free_y')
