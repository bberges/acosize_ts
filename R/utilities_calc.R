computeDistEuclMat <- function(angleMatch,
                           KRMCurrent.current,
                           TS.mat.wide.current,
                           TS.mat.wide.current.info,
                           idx70khz,idx200khz){
  
  ts.model <- KRMCurrent.current$TS
  ts.model <- 10^(ts.model/10)
  
  ts.slow <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'slow',] %>% select(-c('ramping'))
  ts.slow <- 10^(ts.slow/10)
  ts.fast <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'fast',] %>% select(-c('ramping'))
  ts.fast <- 10^(ts.fast/10)
  
  metrics.slow_fast <- as.data.frame(matrix(ncol = 24, nrow = 3))
  colnames(metrics.slow_fast) <- c('angleMatch',
                                   'freqChan',
                                   'minDist',
                                   'sdDist',
                                   'sdDistTrunc25',
                                   'sdDistTrunc10',
                                   'meanDist',
                                   'meanDistTrunc25',
                                   'meanDistTrunc10',
                                   'quant25Dist',
                                   'quant50Dist',
                                   'quant75Dist',
                                   'meanTS.fast',
                                   'meanTS.slow',
                                   'sdTS.fast',
                                   'sdTS.slow',
                                   'meanDistAngleBin.slow',
                                   'meanDistAngleBin.fast',
                                   'minDist.slow',
                                   'minDist.fast',
                                   'idxBestSlow',
                                   'idxBestFast',
                                   'n.fast',
                                   'n.slow')
  
  metrics.model.fast <- as.data.frame(matrix(ncol = 20, nrow = 3))
  colnames(metrics.model.fast) <- c('angleMatch',
                                    'freqChan',
                                    'minDist',
                                    'sdDist',
                                    'sdDistTrunc25',
                                    'sdDistTrunc10',
                                    'meanDist',
                                    'meanDistTrunc25',
                                    'meanDistTrunc10',
                                    'quant25Dist',
                                    'quant50Dist',
                                    'quant75Dist',
                                    'meanTS.mes',
                                    'meanTS.model',
                                    'sdTS.mes',
                                    'meanDistAngleBin',
                                    'minDist.mes',
                                    'idxBest',
                                    'n',
                                    'ramping')
  
  metrics.model.slow <- as.data.frame(matrix(ncol = 20, nrow = 3))
  colnames(metrics.model.slow) <- c('angleMatch',
                                    'freqChan',
                                    'minDist',
                                    'sdDist',
                                    'sdDistTrunc25',
                                    'sdDistTrunc10',
                                    'meanDist',
                                    'meanDistTrunc25',
                                    'meanDistTrunc10',
                                    'quant25Dist',
                                    'quant50Dist',
                                    'quant75Dist',
                                    'meanTS.mes',
                                    'meanTS.model',
                                    'sdTS.mes',
                                    'meanDistAngleBin',
                                    'minDist.mes',
                                    'idxBest',
                                    'n',
                                    'ramping')
  
  matFreqStr <- c('70khz','200khz','all')
  
  for(idx in 1:3){
    if(idx == 1){
      idxSel <- idx70khz
    }else if(idx == 2){
      idxSel <- idx200khz
    }else if(idx == 3){
      idxSel <- idx70khz
      idxSel[1:length(idxSel)] <- TRUE
    }
    
    ts.fast.current <- ts.fast[,idxSel]
    ts.slow.current <- ts.slow[,idxSel]
    ts.model.current <- ts.model[idxSel]
    
    idxFilter.fast <- !apply(is.na(ts.fast.current),1,any)
    ts.fast.current <- ts.fast.current[idxFilter.fast,]
    spectra.str.fast <- subset(TS.mat.wide.current.info,ramping == 'fast')$spectra_id.str[idxFilter.fast]
    
    idxFilter.slow <- !apply(is.na(ts.slow.current),1,any)
    ts.slow.current <- ts.slow.current[idxFilter.slow,]
    spectra.str.slow <- subset(TS.mat.wide.current.info,ramping == 'slow')$spectra_id.str[idxFilter.slow]
    
    # compute metrics for fast ramping
    if(dim(ts.fast.current)[1] != 0){
      dist.fast <- as.vector(dist(ts.fast.current, method = "euclidean"))#/dim(ts.fast.current)[2]
      
      metrics.slow_fast$meanTS.fast[idx]  <- mean(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$sdTS.fast[idx]    <- sd(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$meanDistAngleBin.fast[idx] <- mean(dist.fast,na.rm=TRUE)
      metrics.slow_fast$minDist.fast[idx] <- min(dist.fast,na.rm=TRUE)
    }
    
    # compute metrics for slow ramping
    if(dim(ts.slow.current)[1] != 0){
      dist.slow <- as.vector(dist(ts.slow.current, method = "euclidean"))#/dim(ts.fast.current)[2]
      
      metrics.slow_fast$meanTS.slow[idx]  <- mean(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$sdTS.slow[idx]    <- sd(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$meanDistAngleBin.slow[idx] <- mean(dist.slow,na.rm=TRUE)
      metrics.slow_fast$minDist.slow[idx] <- min(dist.slow,na.rm=TRUE)
    }
    
    # compute metrics between slow and fast ramping
    if(dim(ts.fast.current)[1] != 0 & dim(ts.slow.current)[1] != 0){
      # computing distance matrix metrics between fast and slow ramping
      dist.slow_fast  <- array(0,dim = c(dim(ts.slow.current)[1],dim(ts.fast.current)[1]))
      for(idxSpectra in 1:dim(ts.slow.current)[1]){
        dist.slow_fast[idxSpectra,] <- apply(ts.fast.current,1,TSDistances,y=as.numeric(ts.slow.current[idxSpectra,]),distance="euclidean")
      }
      dist.slow_fast <- dist.slow_fast#/dim(ts.fast.current)[2]
      
      metrics.slow_fast$minDist[idx]          <- min(dist.slow_fast)
      metrics.slow_fast$sdDist[idx]           <- sd(dist.slow_fast)
      metrics.slow_fast$sdDistTrunc25[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$sdDistTrunc10[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$meanDist[idx]             <- mean(dist.slow_fast)
      metrics.slow_fast$meanDistTrunc25[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$meanDistTrunc10[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$quant25Dist[idx]          <- quantile(dist.slow_fast,c(0.25))
      metrics.slow_fast$quant75Dist[idx]          <- quantile(dist.slow_fast,c(0.75))
      metrics.slow_fast$quant50Dist[idx]          <- quantile(dist.slow_fast,c(0.5))
      metrics.slow_fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.slow_fast$idxBestFast[idx]      <- spectra.str.fast[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[2]]
      metrics.slow_fast$idxBestSlow[idx]      <- spectra.str.slow[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[1]]
      metrics.slow_fast$n.fast[idx]           <- dim(ts.fast.current)[1]
      metrics.slow_fast$n.slow[idx]           <- dim(ts.slow.current)[1]
    }
      
    # computing distance matrix metrics between slow/fast ramping and KRM TS model 
    if(dim(ts.fast.current)[1] != 0 & !any(is.na(ts.model.current))){
      dist.fast_model <- apply(ts.fast.current,1,TSDistances,y=ts.model.current,distance="euclidean")
      dist.fast_model <- dist.fast_model#/dim(ts.fast.current)[2]
      
      metrics.model.fast$minDist[idx]          <- min(dist.fast_model)
      metrics.model.fast$sdDist[idx]           <- sd(dist.fast_model)
      metrics.model.fast$sdDistTrunc25[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$sdDistTrunc10[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$meanDist[idx]             <- mean(dist.fast_model)
      metrics.model.fast$meanDistTrunc25[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$meanDistTrunc10[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$quant25Dist[idx]          <- quantile(dist.fast_model,c(0.25))
      metrics.model.fast$quant75Dist[idx]          <- quantile(dist.fast_model,c(0.75))
      metrics.model.fast$quant50Dist[idx]          <- quantile(dist.fast_model,c(0.5))
      metrics.model.fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.fast$idxBest[idx]          <- spectra.str.fast[which(dist.fast_model == min(dist.fast_model), arr.ind = TRUE)]
      metrics.model.fast$ramping[idx]          <- 'fast'
      metrics.model.fast$n[idx]             <- dim(ts.fast.current)[1]
      metrics.model.fast$meanTS.mes[idx]        <- mean(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.fast$meanTS.model[idx]      <- mean(ts.model.current,na.rm=TRUE)
      metrics.model.fast$sdTS.mes[idx]          <- sd(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.fast$meanDistAngleBin[idx]  <- mean(dist.fast,na.rm=TRUE)
      metrics.model.fast$minDist.mes[idx]       <- min(dist.fast,na.rm=TRUE)
    }
    
    if(dim(ts.slow.current)[1] != 0 & !any(is.na(ts.model.current))){
      dist.slow_model <- apply(ts.slow.current,1,TSDistances,y=ts.model.current,distance="euclidean")
      dist.slow_model <- dist.slow_model#/dim(ts.slow.current)[2]
      
      metrics.model.slow$minDist[idx]          <- min(dist.slow_model)
      metrics.model.slow$sdDist[idx]           <- sd(dist.slow_model)
      metrics.model.slow$sdDistTrunc25[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$sdDistTrunc10[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$meanDist[idx]             <- mean(dist.slow_model)
      metrics.model.slow$meanDistTrunc25[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$meanDistTrunc10[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$quant25Dist[idx]          <- quantile(dist.slow_model,c(0.25))
      metrics.model.slow$quant75Dist[idx]          <- quantile(dist.slow_model,c(0.75))
      metrics.model.slow$quant50Dist[idx]          <- quantile(dist.slow_model,c(0.5))
      metrics.model.slow$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.slow$idxBest[idx]          <- spectra.str.slow[which(dist.slow_model == min(dist.slow_model), arr.ind = TRUE)]
      metrics.model.slow$ramping[idx]          <- 'slow'
      metrics.model.slow$n[idx]           <- dim(ts.slow.current)[1]
      metrics.model.slow$meanTS.mes[idx]        <- mean(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.slow$meanTS.model[idx]      <- mean(ts.model.current,na.rm=TRUE)
      metrics.model.slow$sdTS.mes[idx]          <- sd(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.slow$meanDistAngleBin[idx]  <- mean(dist.slow,na.rm=TRUE)
      metrics.model.slow$minDist.mes[idx]       <- min(dist.slow,na.rm=TRUE)
    }
  }
  
  metrics.model <- rbind(metrics.model.slow,metrics.model.fast)
  
  metrics.slow_fast$angleMatch  <- angleMatch
  metrics.model$angleMatch <- angleMatch
  
  return(list(metrics.slow_fast=metrics.slow_fast,
              metrics.model=metrics.model))
}

computeDistDTWMat <- function(angleMatch,
                              KRMCurrent.current,
                              TS.mat.wide.current,
                              TS.mat.wide.current.info,
                              idx70khz,idx200khz){
  ts.model <- KRMCurrent.current$TS
  ts.model <- 10^(ts.model/10)
  
  ts.slow <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'slow',] %>% select(-c('ramping'))
  ts.slow <- 10^(ts.slow/10)
  ts.fast <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'fast',] %>% select(-c('ramping'))
  ts.fast <- 10^(ts.fast/10)
  
  metrics.slow_fast <- as.data.frame(matrix(ncol = 24, nrow = 3))
  colnames(metrics.slow_fast) <- c('angleMatch',
                                   'freqChan',
                                   'minDist',
                                   'sdDist',
                                   'sdDistTrunc25',
                                   'sdDistTrunc10',
                                   'meanDist',
                                   'meanDistTrunc25',
                                   'meanDistTrunc10',
                                   'quant25Dist',
                                   'quant50Dist',
                                   'quant75Dist',
                                   'meanTS.fast',
                                   'meanTS.slow',
                                   'sdTS.fast',
                                   'sdTS.slow',
                                   'meanDistAngleBin.slow',
                                   'meanDistAngleBin.fast',
                                   'minDist.slow',
                                   'minDist.fast',
                                   'idxBestSlow',
                                   'idxBestFast',
                                   'n.fast',
                                   'n.slow')
  
  metrics.model.fast <- as.data.frame(matrix(ncol = 20, nrow = 3))
  colnames(metrics.model.fast) <- c('angleMatch',
                                    'freqChan',
                                    'minDist',
                                    'sdDist',
                                    'sdDistTrunc25',
                                    'sdDistTrunc10',
                                    'meanDist',
                                    'meanDistTrunc25',
                                    'meanDistTrunc10',
                                    'quant25Dist',
                                    'quant50Dist',
                                    'quant75Dist',
                                    'meanTS.mes',
                                    'meanTS.model',
                                    'sdTS.mes',
                                    'meanDistAngleBin',
                                    'minDist.mes',
                                    'idxBest',
                                    'n',
                                    'ramping')
  
  metrics.model.slow <- as.data.frame(matrix(ncol = 20, nrow = 3))
  colnames(metrics.model.slow) <- c('angleMatch',
                                    'freqChan',
                                    'minDist',
                                    'sdDist',
                                    'sdDistTrunc25',
                                    'sdDistTrunc10',
                                    'meanDist',
                                    'meanDistTrunc25',
                                    'meanDistTrunc10',
                                    'quant25Dist',
                                    'quant50Dist',
                                    'quant75Dist',
                                    'meanTS.mes',
                                    'meanTS.model',
                                    'sdTS.mes',
                                    'meanDistAngleBin',
                                    'minDist.mes',
                                    'idxBest',
                                    'n',
                                    'ramping')
  
  matFreqStr <- c('70khz','200khz','all')
  
  for(idx in 1:3){
    if(idx == 1){
      idxSel <- idx70khz
    }else if(idx == 2){
      idxSel <- idx200khz
    }else if(idx == 3){
      idxSel <- idx70khz
      idxSel[1:length(idxSel)] <- TRUE
    }
    
    ts.fast.current <- ts.fast[,idxSel]
    ts.slow.current <- ts.slow[,idxSel]
    ts.model.current <- ts.model[idxSel]
    
    idxFilter.fast <- !apply(is.na(ts.fast.current),1,any)
    ts.fast.current <- ts.fast.current[idxFilter.fast,]
    spectra.str.fast <- subset(TS.mat.wide.current.info,ramping == 'fast')$spectra_id.str[idxFilter.fast]
    
    idxFilter.slow <- !apply(is.na(ts.slow.current),1,any)
    ts.slow.current <- ts.slow.current[idxFilter.slow,]
    spectra.str.slow <- subset(TS.mat.wide.current.info,ramping == 'slow')$spectra_id.str[idxFilter.slow]
    
    # compute metrics for fast ramping
    if(dim(ts.fast.current)[1] != 0){
      list.spectra <- split(as.matrix(ts.fast.current),rep(1:nrow(ts.fast.current),
                                                           each = ncol(ts.fast.current)))
      
      dtw.list <- combn(list.spectra, 2, function(x) {
        dtw(x[[1]], x[[2]],keep.internals = TRUE)}, simplify = FALSE)
      
      dist.fast <- sapply(dtw.list, function(x){as.numeric(x[12])})
      
      metrics.slow_fast$meanTS.fast[idx]  <- mean(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$sdTS.fast[idx]    <- sd(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$meanDistAngleBin.fast[idx] <- mean(dist.fast,na.rm=TRUE)
      metrics.slow_fast$minDist.fast[idx] <- min(dist.fast,na.rm=TRUE)
    }
    
    # compute metrics for slow ramping
    if(dim(ts.slow.current)[1] != 0){
      list.spectra <- split(as.matrix(ts.slow.current),rep(1:nrow(ts.slow.current),
                                                           each = ncol(ts.slow.current)))
      
      dtw.list <- combn(list.spectra, 2, function(x) {
        dtw(x[[1]], x[[2]],keep.internals = TRUE)}, simplify = FALSE)
      
      dist.slow <- sapply(dtw.list, function(x){as.numeric(x[12])})
      
      metrics.slow_fast$meanTS.slow[idx]  <- mean(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$sdTS.slow[idx]    <- sd(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.slow_fast$meanDistAngleBin.slow[idx] <- mean(dist.slow,na.rm=TRUE)
      metrics.slow_fast$minDist.slow[idx] <- min(dist.slow,na.rm=TRUE)
    }
    
    # compute metrics between slow and fast ramping
    if(dim(ts.fast.current)[1] != 0 & dim(ts.slow.current)[1] != 0){
      # computing distance matrix metrics between fast and slow ramping
      dist.slow_fast  <- array(0,dim = c(dim(ts.slow.current)[1],dim(ts.fast.current)[1]))
      for(idxSpectra in 1:dim(ts.slow.current)[1]){
        temp <- apply(ts.fast.current,1,dtw,y=as.numeric(ts.slow.current[idxSpectra,]))
        dist.slow_fast[idxSpectra,] <- sapply(temp, function(x){as.numeric(x[10])})
        #dist.slow_fast[idxSpectra,] <- apply(ts.fast.current,1,TSDistances,y=as.numeric(ts.slow.current[idxSpectra,]),distance="dtw")
      }
      
      # plot(dtw(10*log10(list.spectra[[1]]), 10*log10(list.spectra[[2]]),keep.internals = TRUE),type="two",off=1,match.lty=2,match.indices=20)
      # 
      # windows()
      # plot(dtw(10*log10(as.numeric(ts.slow.current[8,])),
      #          10*log10(as.numeric(ts.fast.current[10,])),
      #          keep=TRUE,step=asymmetric,
      #          open.end=TRUE,open.begin=TRUE),
      #      type="two",off=1,match.lty=2,match.indices=40)
      
      metrics.slow_fast$minDist[idx]          <- min(dist.slow_fast)
      metrics.slow_fast$sdDist[idx]           <- sd(dist.slow_fast)
      metrics.slow_fast$sdDistTrunc25[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$sdDistTrunc10[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$meanDist[idx]             <- mean(dist.slow_fast)
      metrics.slow_fast$meanDistTrunc25[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$meanDistTrunc10[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$quant25Dist[idx]          <- quantile(dist.slow_fast,c(0.25))
      metrics.slow_fast$quant75Dist[idx]          <- quantile(dist.slow_fast,c(0.75))
      metrics.slow_fast$quant50Dist[idx]          <- quantile(dist.slow_fast,c(0.5))
      metrics.slow_fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.slow_fast$idxBestFast[idx]      <- spectra.str.fast[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[2]]
      metrics.slow_fast$idxBestSlow[idx]      <- spectra.str.slow[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[1]]
      metrics.slow_fast$n.fast[idx]           <- dim(ts.fast.current)[1]
      metrics.slow_fast$n.slow[idx]           <- dim(ts.slow.current)[1]
    }
    
    # computing distance matrix metrics between slow/fast ramping and KRM TS model 
    if(dim(ts.fast.current)[1] != 0 & !any(is.na(ts.model.current))){
      temp.dist <- apply(ts.fast.current,1,dtw,y=ts.model.current)
      dist.fast_model <- sapply(temp.dist, function(x){as.numeric(x[10])})
      
      metrics.model.fast$minDist[idx]          <- min(dist.fast_model)
      metrics.model.fast$sdDist[idx]           <- sd(dist.fast_model)
      metrics.model.fast$sdDistTrunc25[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$sdDistTrunc10[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$meanDist[idx]             <- mean(dist.fast_model)
      metrics.model.fast$meanDistTrunc25[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$meanDistTrunc10[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$quant25Dist[idx]          <- quantile(dist.fast_model,c(0.25))
      metrics.model.fast$quant75Dist[idx]          <- quantile(dist.fast_model,c(0.75))
      metrics.model.fast$quant50Dist[idx]          <- quantile(dist.fast_model,c(0.5))
      metrics.model.fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.fast$idxBest[idx]          <- spectra.str.fast[which(dist.fast_model == min(dist.fast_model), arr.ind = TRUE)]
      metrics.model.fast$ramping[idx]          <- 'fast'
      metrics.model.fast$n[idx]             <- dim(ts.fast.current)[1]
      metrics.model.fast$meanTS.mes[idx]        <- mean(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.fast$meanTS.model[idx]      <- mean(ts.model.current,na.rm=TRUE)
      metrics.model.fast$sdTS.mes[idx]          <- sd(apply(ts.fast.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.fast$meanDistAngleBin[idx]  <- mean(dist.fast,na.rm=TRUE)
      metrics.model.fast$minDist.mes[idx]       <- min(dist.fast,na.rm=TRUE)
    }
    
    if(dim(ts.slow.current)[1] != 0 & !any(is.na(ts.model.current))){
      temp.dist <- apply(ts.slow.current,1,dtw,y=ts.model.current)
      dist.slow_model <- sapply(temp.dist, function(x){as.numeric(x[10])})
      
      metrics.model.slow$minDist[idx]          <- min(dist.slow_model)
      metrics.model.slow$sdDist[idx]           <- sd(dist.slow_model)
      metrics.model.slow$sdDistTrunc25[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$sdDistTrunc10[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$meanDist[idx]             <- mean(dist.slow_model)
      metrics.model.slow$meanDistTrunc25[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$meanDistTrunc10[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$quant25Dist[idx]          <- quantile(dist.slow_model,c(0.25))
      metrics.model.slow$quant75Dist[idx]          <- quantile(dist.slow_model,c(0.75))
      metrics.model.slow$quant50Dist[idx]          <- quantile(dist.slow_model,c(0.5))
      metrics.model.slow$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.slow$idxBest[idx]          <- spectra.str.slow[which(dist.slow_model == min(dist.slow_model), arr.ind = TRUE)]
      metrics.model.slow$ramping[idx]          <- 'slow'
      metrics.model.slow$n[idx]           <- dim(ts.slow.current)[1]
      metrics.model.slow$meanTS.mes[idx]        <- mean(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.slow$meanTS.model[idx]      <- mean(ts.model.current,na.rm=TRUE)
      metrics.model.slow$sdTS.mes[idx]          <- sd(apply(ts.slow.current,1,mean,na.rm=TRUE),na.rm=TRUE)
      metrics.model.slow$meanDistAngleBin[idx]  <- mean(dist.slow,na.rm=TRUE)
      metrics.model.slow$minDist.mes[idx]       <- min(dist.slow,na.rm=TRUE)
    }
  }
  
  metrics.model <- rbind(metrics.model.slow,metrics.model.fast)
  
  metrics.slow_fast$angleMatch  <- angleMatch
  metrics.model$angleMatch <- angleMatch
  
  return(list(metrics.slow_fast=metrics.slow_fast,
              metrics.model=metrics.model))
}

computeDistDTWMat.old <- function(angleMatch,
                               KRMCurrent.current,
                               TS.mat.wide.current,
                               TS.mat.wide.current.info,
                               idx70khz,idx200khz){
  ts.model <- KRMCurrent.current$TS
  
  ts.slow <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'slow',] %>% select(-c('ramping'))
  ts.fast <- TS.mat.wide.current[TS.mat.wide.current$ramping == 'fast',] %>% select(-c('ramping'))
  
  metrics.slow_fast <- as.data.frame(matrix(ncol = 16, nrow = 3))
  colnames(metrics.slow_fast) <- c('angleMatch','freqChan','minDist','sdDist','sdDistTrunc25','sdDistTrunc10','mean','meanDistTrunc25','meanDistTrunc10','quant25','quant50','quant75','idxBestSlow','idxBestFast','n.fast','n.slow')
  
  metrics.model.fast <- as.data.frame(matrix(ncol = 15, nrow = 3))
  colnames(metrics.model.fast) <- c('angleMatch','freqChan','minDist','sdDist','sdDistTrunc25','sdDistTrunc10','mean','meanDistTrunc25','meanDistTrunc10','quant25','quant50','quant75','idxBest','n','ramping')
  
  metrics.model.slow <- as.data.frame(matrix(ncol = 15, nrow = 3))
  colnames(metrics.model.slow) <- c('angleMatch','freqChan','minDist','sdDist','sdDistTrunc25','sdDistTrunc10','mean','meanDistTrunc25','meanDistTrunc10','quant25','quant50','quant75','idxBest','n','ramping')
  
  matFreqStr <- c('70khz','200khz','all')
  
  for(idx in 1:3){
    if(idx == 1){
      idxSel <- idx70khz
    }else if(idx == 2){
      idxSel <- c(idx200khz)
    }else if(idx == 3){
      idxSel <- idx70khz
      idxSel[1:length(idxSel)] <- TRUE
    }
    
    ts.fast.current <- ts.fast[,idxSel]
    ts.slow.current <- ts.slow[,idxSel]
    ts.model.current <- ts.model[idxSel]
    
    idxFilter.fast <- !apply(is.na(ts.fast.current),1,any)
    ts.fast.current <- ts.fast.current[idxFilter.fast,]
    spectra.str.fast <- subset(TS.mat.wide.current.info,ramping == 'fast')$spectra_id.str[idxFilter.fast]
    
    idxFilter.slow <- !apply(is.na(ts.slow.current),1,any)
    ts.slow.current <- ts.slow.current[idxFilter.slow,]
    spectra.str.slow <- subset(TS.mat.wide.current.info,ramping == 'slow')$spectra_id.str[idxFilter.slow]
    
    if(dim(ts.fast.current)[1] != 0 & dim(ts.slow.current)[1] != 0){
      # computing distance matrix metrics between fast and slow ramping
      dist.slow_fast  <- array(0,dim = c(dim(ts.slow.current)[1],dim(ts.fast.current)[1]))
      for(idxSpectra in 1:dim(ts.slow.current)[1]){
        temp <- apply(ts.fast.current,1,dtw,y=as.numeric(ts.slow.current[idxSpectra,]))
        dist.slow_fast[idxSpectra,] <- sapply(temp, function(x){as.numeric(x[10])})
        #dist.slow_fast[idxSpectra,] <- apply(ts.fast.current,1,TSDistances,y=as.numeric(ts.slow.current[idxSpectra,]),distance="dtw")
      }
      
      output <- dtw(as.numeric(ts.fast.current[1,]),as.numeric(ts.slow.current[idxSpectra,]))
      
      metrics.slow_fast$minDist[idx]          <- min(dist.slow_fast)
      metrics.slow_fast$sdDist[idx]           <- sd(dist.slow_fast)
      metrics.slow_fast$sdDistTrunc25[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$sdDistTrunc10[idx]    <- sd(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$mean[idx]             <- mean(dist.slow_fast)
      metrics.slow_fast$meanDistTrunc25[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.75))])
      metrics.slow_fast$meanDistTrunc10[idx]  <- mean(dist.slow_fast[dist.slow_fast < quantile(dist.slow_fast,c(0.90))])
      metrics.slow_fast$quant25[idx]          <- quantile(dist.slow_fast,c(0.25))
      metrics.slow_fast$quant75[idx]          <- quantile(dist.slow_fast,c(0.75))
      metrics.slow_fast$quant50[idx]          <- quantile(dist.slow_fast,c(0.5))
      metrics.slow_fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.slow_fast$idxBestFast[idx]      <- spectra.str.fast[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[2]]
      metrics.slow_fast$idxBestSlow[idx]      <- spectra.str.slow[which(dist.slow_fast == min(dist.slow_fast), arr.ind = TRUE)[1]]
      metrics.slow_fast$n.fast[idx]           <- dim(ts.fast.current)[1]
      metrics.slow_fast$n.slow[idx]           <- dim(ts.slow.current)[1]
    }
    
    # computing distance matrix metrics between slow/fast ramping and KRM TS model 
    if(dim(ts.fast.current)[1] != 0 & !any(is.na(ts.model.current))){
      temp <- apply(ts.fast.current,1,dtw,y=ts.model.current)
      dist.fast_model <- sapply(temp, function(x){as.numeric(x[10])})
      
      #dist.fast_model <- apply(ts.fast.current,1,TSDistances,y=ts.model.current,distance="dtw")
      
      metrics.model.fast$minDist[idx]          <- min(dist.fast_model)
      metrics.model.fast$sdDist[idx]           <- sd(dist.fast_model)
      metrics.model.fast$sdDistTrunc25[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$sdDistTrunc10[idx]    <- sd(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$mean[idx]             <- mean(dist.fast_model)
      metrics.model.fast$meanDistTrunc25[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.75))])
      metrics.model.fast$meanDistTrunc10[idx]  <- mean(dist.fast_model[dist.fast_model < quantile(dist.fast_model,c(0.90))])
      metrics.model.fast$quant25[idx]          <- quantile(dist.fast_model,c(0.25))
      metrics.model.fast$quant75[idx]          <- quantile(dist.fast_model,c(0.75))
      metrics.model.fast$quant50[idx]          <- quantile(dist.fast_model,c(0.5))
      metrics.model.fast$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.fast$idxBest[idx]          <- spectra.str.fast[which(dist.fast_model == min(dist.fast_model), arr.ind = TRUE)]
      metrics.model.fast$ramping[idx]          <- 'fast'
      metrics.model.fast$n[idx]           <- dim(ts.fast.current)[1]
    }
    
    if(dim(ts.slow.current)[1] != 0 & !any(is.na(ts.model.current))){
      temp <- apply(ts.slow.current,1,dtw,y=ts.model.current)
      dist.slow_model <- sapply(temp, function(x){as.numeric(x[10])})
      
      #dist.slow_model <- apply(ts.slow.current,1,TSDistances,y=ts.model.current,distance="dtw")
      
      metrics.model.slow$minDist[idx]          <- min(dist.slow_model)
      metrics.model.slow$sdDist[idx]           <- sd(dist.slow_model)
      metrics.model.slow$sdDistTrunc25[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$sdDistTrunc10[idx]    <- sd(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$mean[idx]             <- mean(dist.slow_model)
      metrics.model.slow$meanDistTrunc25[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.75))])
      metrics.model.slow$meanDistTrunc10[idx]  <- mean(dist.slow_model[dist.slow_model < quantile(dist.slow_model,c(0.90))])
      metrics.model.slow$quant25[idx]          <- quantile(dist.slow_model,c(0.25))
      metrics.model.slow$quant75[idx]          <- quantile(dist.slow_model,c(0.75))
      metrics.model.slow$quant50[idx]          <- quantile(dist.slow_model,c(0.5))
      metrics.model.slow$freqChan[idx]         <- matFreqStr[idx]
      metrics.model.slow$idxBest[idx]          <- spectra.str.slow[which(dist.slow_model == min(dist.slow_model), arr.ind = TRUE)]
      metrics.model.slow$ramping[idx]          <- 'slow'
      metrics.model.slow$n[idx]           <- dim(ts.slow.current)[1]
    }
  }
  
  metrics.model <- rbind(metrics.model.slow,metrics.model.fast)
  
  metrics.slow_fast$angleMatch  <- angleMatch
  metrics.model$angleMatch <- angleMatch
  
  return(list(metrics.slow_fast=metrics.slow_fast,
              metrics.model=metrics.model))
}

#Get volume
#V = pi * a * b * h
getV = function(x,zL,zU,w){
  sum(na.omit(c(0,diff(x)) * (zU-zL) * w ))
}