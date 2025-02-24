#path <- 'J:/git/heras_index_kbwot'
path <- 'C:/git/acosize_ts'

try(setwd(path),silent=TRUE)

myPath <- data.frame(mainPath=file.path("."),
                     dataPath=file.path(".","data"),
                     figuresPath=file.path(".","figures"),
                     resultsPath=file.path(".","results"),
                     softwarePath=file.path(".","software"))

#install.packages(file.path(myPath$softwarePath,'KRMr_0.4.2.tar.gz'), repos = NULL, type="source")

library(KRMr)
library(imager)
library(ggplot2)

lengthData <- read.csv(file.path(myPath$dataPath,'length_data.csv'))
summaryTab <- read.csv(file.path(myPath$resultsPath,'summaryTab.csv'))
summaryTab <- summaryTab[summaryTab$KRM == 1,]
summaryTab <- summaryTab[summaryTab$type != 'MAC',]

#lengthData <- lengthData[lengthData$fish_id %in% unique(summaryTab$fish_id),]

#shapes.mac <- readRDS(file.path(myPath$dataPath,'xrays','shapes','shapes.mac.RDS'))

frequency=(12:300)*1e3 #frequencies Hz
theta = c(seq(65,115)) #orientations

c.w = 1490 #Ambient water m/s
rho.w = 1026.8 #density ambient water kg/m3

c.fb=1570 #soundspeed body 
c.b=345 #soundspeed swimbladder
cs = c(c.fb, c.b)

rho.fb=1070 #density body kg/m3
rho.b=1.24 #density swimbladder kg/m3
rhos = c(rho.fb, rho.b)

rotAngle <- pi

# lh.rel = c(0.93493081,1.085137085,0.986937831,0.992118607,1.043363239,0.971804511,
#            1.008151008,0.943510738,0.888715413,0.962272832,0.998366013,0.971395503,
#            0.859066859,1.091848545,0.876022126,0.952594299,0.903527337,1,1.064192344,
#            0.881422693,0.962962963,0.883969907,1.038593215,0.90206453) #Relative scaling based on infor in the length Excel sheet
# dy = shape$z_U - shape$z_L
# dym= (shape$z_U + shape$z_L)/2
# shapes.mac= lapply(lh.rel, FUN=function(sc){
#   data.frame(x=shape$x,
#              w=shape$w,
#              z_U=dym+sc*(dy/2),
#              z_L=dym-sc*(dy/2))})
# 
# 
# shp01=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0('M10','.csv')))
# shp02=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0('M12','.csv')))
# 
# shp01_lateral <- subset(shp01,asp1 == 'Lateral_body')
# max(shp01_lateral$BX)-min(shp01_lateral$BX)
# 
# shp02_lateral <- subset(shp02,asp1 == 'Lateral_body')
# max(shp02_lateral$BX)-min(shp02_lateral$BX)

for(fish_id in unique(summaryTab$fish_id)){
  print(fish_id)
  currentL    <- lengthData$length_2[lengthData$fish_id == fish_id]*1e-2
  shapeMatch  <- unique(summaryTab[summaryTab$fish_id == fish_id,]$shapeMatch)
  
  if(fish_id == 'S23' & substr(fish_id,1,1) != 'M'){
    shp=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0(fish_id,'.csv')))
    lb=shp[shp$asp1=="Lateral_bladder",]
    shp=as.data.frame(shp)%>%dplyr::filter(asp1!="Lateral_bladder")%>%rbind(lb[c(22:99,seq(21,1,by=-1)),])
    currentShpRot <- Imagej2shp(shp=shp,
                            dorsal=c("Lateral_body","Lateral_bladder"), 
                            lateral = c("Dorsal_body","Dorsal_bladder"))
    # currentShp = Imagej2shp(shp=shp,  dorsal=c("Dorsal_body","Dorsal_bladder"), lateral=c("Lateral_body", "Lateral_bladder"))
    currentShp <- currentShpRot
    for(idx in 1:length(currentShpRot)){
      currentShp[[idx]]$x    <- currentShp[[idx]]$x*cos(rotAngle)#+shpRot$BY*sin(rotAngle)
      currentShp[[idx]]$z_U  <- -currentShp[[idx]]$x*sin(rotAngle)-currentShp[[idx]]$z_U*cos(rotAngle)
      currentShp[[idx]]$z_L  <- -currentShp[[idx]]$x*sin(rotAngle)-currentShp[[idx]]$z_L*cos(rotAngle)
    }
  }else if(substr(fish_id,1,1) != 'M'){
    shp=read.csv(file.path(myPath$dataPath,'xrays','shapes',paste0(fish_id,'.csv')))
    currentShpRot = Imagej2shp(shp=shp,
                               dorsal=c("Lateral_body","Lateral_bladder"), 
                               lateral = c("Dorsal_body","Dorsal_bladder"))
    currentShp <- currentShpRot
    for(idx in 1:length(currentShpRot)){
      currentShp[[idx]]$x    <- currentShp[[idx]]$x*cos(rotAngle)#+shpRot$BY*sin(rotAngle)
      currentShp[[idx]]$z_U  <- -currentShp[[idx]]$x*sin(rotAngle)-currentShp[[idx]]$z_U*cos(rotAngle)
      currentShp[[idx]]$z_L  <- -currentShp[[idx]]$x*sin(rotAngle)-currentShp[[idx]]$z_L*cos(rotAngle)
    }
  }else if(substr(fish_id,1,1) == 'M'){
    idxMAC <- as.numeric(strsplit(fish_id,'M')[[1]][2])
    currentShp <- shapes.mac[[idxMAC]]
    currentShpRot <- currentShp
    currentShpRot$x    <- currentShpRot$x*cos(rotAngle)#+shpRot$BY*sin(rotAngle)
    currentShpRot$z_U  <- -currentShpRot$x*sin(rotAngle)-currentShpRot$z_U*cos(rotAngle)
    currentShpRot$z_L  <- -currentShpRot$x*sin(rotAngle)-currentShpRot$z_L*cos(rotAngle)
  }
  
  # M10L <- max(subset(shpM10,asp1 == 'Lateral_body')$BX)-min(subset(shpM10,asp1 == 'Lateral_body')$BX)
  # 
  # currentL/(M10L*1e-3)
  # 
  # dy = shape$z_U - shape$z_L
  # dym= (shape$z_U + shape$z_L)/2
  # shapes.mac= lapply(lh.rel, FUN=function(sc){
  #   data.frame(x=shape$x,
  #              w=shape$w,
  #              z_U=dym+sc*(dy/2),
  #              z_L=dym-sc*(dy/2))})
  
  
  # KRMr::get_shp3d(currentShp)
  # KRMr::get_shp3d(currentShpRot)
  
  # ggplot()+
  #   geom_line(data=currentShp[[1]],aes(x=x,y=z_U))+
  #   geom_line(data=currentShp[[1]],aes(x=x,y=z_L))+
  #   geom_line(data=currentShp[[2]],aes(x=x,y=z_U))+
  #   geom_line(data=currentShp[[2]],aes(x=x,y=z_L))+
  #   geom_line(data=currentShpRot[[1]],aes(x=x,y=z_U,col='rotated'))+
  #   geom_line(data=currentShpRot[[1]],aes(x=x,y=z_L,col='rotated'))+
  #   geom_line(data=currentShpRot[[2]],aes(x=x,y=z_U,col='rotated'))+
  #   geom_line(data=currentShpRot[[2]],aes(x=x,y=z_L,col='rotated'))

  if(substr(fish_id,1,1) == 'M'){
    currentKRM = krm(frequency = frequency,
                     c.w = c.w,
                     rho.w = rho.w,
                     theta = theta,
                     cs =1570,
                     rhos = 1060,
                     shape = currentShp,
                     L = currentL,
                     fb = 1)
    
    currentKRMRot = krm(frequency = frequency,
                        c.w = c.w,
                        rho.w = rho.w,
                        theta = theta,
                        cs =1570,
                        rhos = 1060,
                        shape = currentShpRot,
                        L = currentL,
                        fb = 1)
  }else if(substr(fish_id,1,1) != 'M'){
    currentKRM = krm(frequency = frequency,
                     c.w = c.w,
                     rho.w = rho.w,
                     theta = theta,
                     cs =cs,
                     rhos = rhos,
                     shape = currentShp,
                     modes = c("fluid","soft"),
                     L = currentL)
    
    currentKRMRot = krm(frequency = frequency,
                        c.w = c.w,
                        rho.w = rho.w,
                        theta = theta,
                        cs =cs,
                        rhos = rhos,
                        shape = currentShpRot,
                        modes = c("fluid","soft"),
                        L = currentL,
                        fb = 1)
  }
  
  currentKRMRot$theta <- currentKRMRot$theta+180
  currentKRM <- rbind(currentKRMRot,currentKRM)
  
  # windows()
  # ggplot(subset(currentKRM,theta == 115 | theta == 295),aes(x=frequency,y=TS,col=as.factor(theta)))+
  #   geom_line()
  
  write.csv(currentKRM,file = file.path(myPath$resultsPath,'KRM',paste0('krm_',fish_id,'.csv')),row.names = FALSE)
}
