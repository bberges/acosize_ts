#path <- 'J:/git/heras_index_kbwot'
path <- 'G:/git/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     softwarePath=file.path(".","software"))

install.packages(file.path(myPath$softwarePath,'KRMr_0.4.2.tar.gz'), repos = NULL, type="source")

library(KRMr)
library(imager)
library(ggplot2)

lengthData <- read.csv(file.path(myPath$dataPath,'length_data.csv'))
summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]

lengthData <- lengthData[lengthData$fish_id %in% unique(summaryTab$fish_id),]

frequency=(12:300)*1e3 #frequencies Hz
theta = seq(65,115) #orientations

c.w = 1490 #Ambient water m/s
rho.w = 1026.8 #density ambient water kg/m3

c.fb=1570 #soundspeed body 
c.b=345 #soundspeed swimbladder
cs = c(c.fb, c.b)

rho.fb=1070 #density body kg/m3
rho.b=1.24 #density swimbladder kg/m3
rhos = c(rho.fb, rho.b)



for(fish_id in unique(summaryTab$fish_id)[28:31]){
  print(fish_id)
  currentL <- lengthData$length_2[lengthData$fish_id == fish_id]
  currentShp = Imagej2shp(shp=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0(fish_id,'.csv'))))
  
  currentKRM = krm(frequency = frequency,
                   c.w = c.w,
                   rho.w = rho.w,
                   theta = theta,
                   cs =cs,
                   rhos = rhos,
                   shape = currentShp,
                   modes = c("fluid","soft"),
                   L = currentL,
                   fb = 1)
  
  write.csv(currentKRM,file = file.path(myPath$resultsPath,'KRM',paste0('krm_',fish_id,'.csv')))
}
