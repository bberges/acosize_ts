rm(list = ls())

#### Packages ####
library(KRMr)
library(imager)
library(ggplot2)
library(tidyverse)
library(icesTAF)
library(doParallel)


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

ncores <- detectCores()-1
cl <- makeCluster(ncores)
registerDoParallel(cl)
clusterEvalQ(cl,library(KRMr))
clusterEvalQ(cl,library(ggplot2))
clusterEvalQ(cl,library(tidyverse))
clusterEvalQ(cl,source('G:/git/acosize_ts/R/utilities_plot.R'))

foreach(idxAnimal = uniqueAnimals) %dopar% plotFun(myPath,summaryTab,idxAnimal)

