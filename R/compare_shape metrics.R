#path <- 'J:/git/heras_index_kbwot'
path <- 'C:/git/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     softwarePath=file.path(".","software"),
                     rPath=file.path(".","R"))

#install.packages(file.path(myPath$softwarePath,'KRMr_0.4.2.tar.gz'), repos = NULL, type="source")

library(KRMr)
library(imager)
library(ggplot2)

source(file.path(myPath$rPath,'utilities_calc.R'))

lengthData <- read.csv(file.path(myPath$dataPath,'length_data.csv'))
summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]

#lengthData <- lengthData[lengthData$fish_id %in% unique(summaryTab$fish_id),]

shapes.mac <- readRDS(file.path(myPath$dataPath,'xrays','shapes','shapes.mac.RDS'))

uniqueAnimals <- unique(summaryTab$fish_id)

metrics.shape <- as.data.frame(matrix(ncol = 9, nrow = length(uniqueAnimals)))
colnames(metrics.shape) <- c('fish_id',
                             'length',
                             'width',
                             'bodyV',
                             'bladderV',
                             'ratioV',
                             'rationWL',
                             'species',
                             'ratioLV')
metrics.shape$fish_id <- uniqueAnimals

for(fish_id in uniqueAnimals){
  print(fish_id)
  idxCurrent <- which(metrics.shape$fish_id == fish_id)
  currentL <- lengthData$length_2[lengthData$fish_id == fish_id]*1e-2
  currentW <- lengthData$width[lengthData$fish_id == fish_id]*1e-2
  species <- lengthData$species[lengthData$fish_id == fish_id]
  
  if(fish_id == 'S23' & substr(fish_id,1,1) != 'M'){
    shp=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0(fish_id,'.csv')))
    lb=shp[shp$asp1=="Lateral_bladder",]
    shp=as.data.frame(shp)%>%dplyr::filter(asp1!="Lateral_bladder")%>%rbind(lb[c(22:99,seq(21,1,by=-1)),])
    currentShp = Imagej2shp(shp=shp,  dorsal=c("Dorsal_body","Dorsal_bladder"), lateral=c("Lateral_body", "Lateral_bladder"))
  }else if(substr(fish_id,1,1) != 'M'){
    shp=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0(fish_id,'.csv')))
    currentShp <- Imagej2shp(shp=shp)
  }else if(substr(fish_id,1,1) == 'M'){
    idxMAC <- as.numeric(strsplit(fish_id,'M')[[1]][2])
    currentShp <- shapes.mac[[idxMAC]]
  }
  
  if(substr(fish_id,1,1) != 'M'){
    bladderV <- getV(currentShp[[2]]$x,
                     currentShp[[2]]$z_L,
                     currentShp[[2]]$z_U,
                     currentShp[[2]]$w)
    
    bodyV <- getV(currentShp[[1]]$x,
                  currentShp[[1]]$z_L,
                  currentShp[[1]]$z_U,
                  currentShp[[1]]$w)
    
    metrics.shape$bladderV[idxCurrent]  <- bladderV
    metrics.shape$ratioV[idxCurrent]    <- bladderV/bodyV
  }else{
    bodyV <- getV(currentShp$x,
                  currentShp$z_L,
                  currentShp$z_U,
                  currentShp$w)
  }
  
  metrics.shape$length[idxCurrent]    <- currentL
  metrics.shape$width[idxCurrent]     <- currentW
  metrics.shape$bodyV[idxCurrent]     <- bodyV
  metrics.shape$ratioLV[idxCurrent]   <- currentL/bodyV
  metrics.shape$rationWL[idxCurrent]  <- currentW/currentL
  metrics.shape$species[idxCurrent]   <- species
}

save(metrics.shape,
     file = file.path(myPath$resultsPath,paste0('fish_dimensions.RData')))

###############################################################################
# plotting
###############################################################################

scaling_factor <- 2
png(file.path(myPath$figuresPath,'fish_length_histogram.png'), 
    width = 8*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(metrics.shape,aes(x=length*1e2,fill=species))+
      geom_histogram(binwidth=1)+
      xlab('Fish length (cm)')+
      ylab('Count (#)')

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'fish_ratioV.png'), 
    width = 8*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.shape,!is.na(ratioV)),aes(y=ratioV*100,x=species))+
      geom_boxplot()+
      geom_jitter(aes(col=species))+
  ylab('Swimbladder/body ratio (%)')+
  theme(legend.position = 'none')

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'fish_ratioWL.png'), 
    width = 8*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(subset(metrics.shape,!is.na(rationWL)),aes(y=rationWL*100,x=species))+
      geom_boxplot()+
      geom_jitter(aes(col=species))+
      ylab('width/length (%)')+
      theme(legend.position = 'none')

print(p)
dev.off()

scaling_factor <- 2
png(file.path(myPath$figuresPath,'fish_width.png'), 
    width = 8*scaling_factor, height = 8*scaling_factor, units = "cm", res = 300, pointsize = 10)

p <- ggplot(metrics.shape,aes(y=width,x=species))+
  geom_boxplot()+
  geom_jitter(aes(col=species))+
  ylab('Width (cm)')+
  theme(legend.position = 'none')

print(p)
dev.off()