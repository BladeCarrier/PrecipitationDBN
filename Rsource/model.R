

library("caret")
library("deepnet")
library(foreach)
library(doSNOW)

source("Rsource/model.base.R")
source("Rsource/model.auxilary.R")
source("Rsource/model.serial.R")




model.getLinearModels.arlm = function(model.cfg,layer.cfg, rbm.cfg, RBMups, andsave = FALSE)
{

  
  reg.L2 = layer.cfg$reg.L2
  order = layer.cfg$order
  
  numThreads = model.cfg$numThreads
  rbm.cfg = model.cfg$rbm.l0
  numHidden = rbm.cfg$numHidden
  sectorWidth = rbm.cfg$sector.width
  sectorHeight = rbm.cfg$sector.height

  toXData = rbm.toXData
  toYData = rbm.toYData

  modelLonDiap = rbm.cfg$layerLonDiap
  modelLatDiap = rbm.cfg$layerLatDiap
  
  linearModelMatrices = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(linearModelMatrices) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
  encoder = rbm.cfg$encoder
  encode = encoder.encode
  
  logit = linearModel.logit

  
  gc()
  
  cl <- makeCluster(numThreads)
  registerDoSNOW(cl)
  print("training LMs")
  for (lon in modelLonDiap)
  {
    print(paste("now processing lon", lon))
    lres = foreach(lat = modelLatDiap,.combine = "cbind") %dopar%
    {
      library(caret)
      library(foreach)
      library(penalized)
      library(stepPlr)
      
      up = RBMups[[paste("tile.",lon,'.',lat,sep = "")]]
      X=(order+1):dim(RBMups)[1]
    
      xdata = foreach(day = X, .combine = "rbind") %do% toXData(numHidden,up,order,day)
      xdata = as.data.frame.array(xdata) 
      
      ydata = foreach(day = X, .combine = "rbind") %do% toYData(numHidden,up,day)
      ydata = as.data.frame.array(ydata)
      
      cMatr = foreach(yimg = 1:numHidden, .combine = "rbind") %do%
      {
        #dframe = xdata
        #dframe$targ = as.numeric( ydata[,paste("y",yimg)])
        
        cVec = tryCatch({
          label = as.numeric( ydata[,paste("y",yimg)])
          linear =plr(y = label, x = xdata,lambda=reg.L2)
          coefficients( linear,"all")
          
        }, error = function(w) {
          print(w)
          interc = logit(mean(dframe$targ))
          as.vector(c(interc,rep(0,ncol(xdata))))
        })
        
        
        
        cVec
      }

      remove(xdata)
      remove(ydata)
      gc()
      c(paste("tile.",lon,'.',lat,sep = ""), list(cMatr))
    }
    for (r in 1:ncol(lres))
    {
      linearModelMatrices[[lres[[1,r]]]] = as.matrix(lres[[2,r]])
    }
    rm(lres)
  }

  stopCluster(cl)
  if (andsave)
    model.save(linearModelMatrices,model.cfg,".arlms")
  return(linearModelMatrices)
}


model.makePredictions.arlm = function(model.cfg,layer.cfg, rbm.cfg,rbmup,arlms,dims)
{
  #predicts for every day between [order+1]th and [days in RBMup+1]th
  print("predicting with linear models")
  require(foreach)
  require(doSNOW)
  require(deepnet)
  

  order = layer.cfg$order
  
  rbm.cfg = model.cfg$rbm.l0
  numHidden = rbm.cfg$numHidden
  sectorWidth = rbm.cfg$sector.width
  sectorHeight = rbm.cfg$sector.height
  
    
  toXData = rbm.toXData
  lmpredict = linearModel.predict
  getDiap = model.grid.getDiap
  
  modelLonDiap = rbm.cfg$layerLonDiap
  modelLatDiap = rbm.cfg$layerLatDiap
  
  RBMpreds = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(RBMpreds) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
  
  dayDiap = (order+1):(dim(rbmup)[1]+1) #days FOR WHICH a prediction is made
  
  batchSize = model.cfg$batchSize
  numThreads = model.cfg$numThreads
  
  pairCoords = expand.grid(modelLonDiap,modelLatDiap)
  pairCoords = pairCoords[order(pairCoords[,1]), ]
  
  cl <- makeCluster(numThreads)
  registerDoSNOW(cl)
  currLon = 0#is only used for progress reports
  while(TRUE)
  {
    batch = pairCoords[1:batchSize,]
    batch = batch[!is.na(batch[,1]),]
    
    if (batch[1,1] != currLon)
    {#report status
      print(paste("processing lon",batch[1,1]))
      currLon = batch[1,1]
    }
    
    rbmUpLocal = list()
    arlmLocal = list()
    for(i in 1:dim(batch)[1]) 
    {
      rbmind = paste("tile.",batch[i,1],'.',batch[i,2],sep="")
      rbmUpLocal[[rbmind]] = rbmup[[rbmind]] 
      arlmLocal[[rbmind]] = arlms[[rbmind]]
    }
    
    
    lres = foreach(i = 1:dim(batch)[1],.combine = "cbind",.noexport = "arlms") %dopar%
    {
      library(deepnet)
      library(foreach)
      lon = batch[i,1]
      lat = batch[i,2]
      #check if this works for several days of rbmup at once
      lonDiap = getDiap(lon,sectorWidth,dims[2])
      latDiap = getDiap(lat,sectorHeight,dims[3])
      
      arlm = arlmLocal[[paste("tile.",lon,'.',lat,sep="")]]
      up = rbmUpLocal[[paste("tile.",lon,'.',lat,sep="")]]
      
      predRes = foreach(day = dayDiap,.combine="rbind") %do%
      {
        lastUp = up[(day - order):(day-1),]
        datarow = toXData(numHidden,lastUp,order,day =  order+1)#predict for next day
        output = lmpredict(arlm,datarow)#RBMDOWN-ready?
        output
      
      }
      gc()  
      c(paste("tile.",lon,'.',lat,sep=""), list(predRes))
      

    }

    dim(lres) = c(2,nrow(batch))#eliminating one-element batch case (ncol = NULL)
    for (r in 1:ncol(lres))
    {
      RBMpreds[[lres[[1,r]]]] = matrix(lres[[2,r]],nrow = length(dayDiap),ncol = length(lres[[2,r]])/length(dayDiap))
    }
    
    
    if (batchSize >= dim(pairCoords)[1])
      break
    pairCoords = pairCoords[(batchSize+1):dim(pairCoords)[1],]

  }
  stopCluster(cl)

  gc()
  return(RBMpreds)
}
model.predict = function(model, pratemap)
{
  dims = dim(pratemap)
  pratemap.matr = matrix(pratemap,nrow = dims[1], ncol = dims[2]*dims[3])
  pratemap.diff = diff(pratemap.matr)
  dim(pratemap.diff) = dim(pratemap) - c(1,0,0)
  rm(pratemap.matr)
  
  RBMups = model.makeRBMup.base(model$rbm.l0,
                                model$cfg,
                                model$cfg$rbm.l0,
                                pratemap.diff)
  RBMpreds = model.makePredictions.arlm(model$cfg,
                                        model$cfg$pred.l0,
                                        model$cfg$rbm.l0,
                                        RBMups,
                                        model$arlms,
                                        dim(pratemap.diff))
  
  dims = dim(pratemap) - c(model$cfg$pred.l0$order,0,0)
  RBMdowns = model.makeRBMdown.base(model$rbm.l0,
                                    model$cfg,
                                    model$cfg$rbm.l0,
                                    RBMpreds,
                                    dims  )
  #RBMdowns are predictions of difference between ([order+1] to [order+2]) and ([numDays] to [numDays+1])
  order = model$cfg$pred.l0$order
  dims = dim(pratemap) + c(1,0,0)
  ans = array(NA, dim = dims)
  
  #prediction for ith day is sum of the known previous day plus predicted difference between previous and ith.
  ans[(order+2):dims[1],,] = pratemap[(order+1):dim(pratemap)[1],,] +RBMdowns
  print("note that several first days will be NA because model needs a certain amount of history to predict")
  return(ans)
  
}
model.save = function(object,cfg,mark = "")
{
  print("saving...")
  cVec = tryCatch({
    modelFileName = paste(cfg$path,cfg$name,".RData",sep = "")
    save(object,file = modelFileName)
    print(paste("model saved to",modelFileName))
    
  }, error = function(w) {
    print(w)
    print("failed to save model")
  })
}
model.assemble = function(cfg,pratemap,andsave = FALSE)
{
  dims = dim(pratemap)
  pratemap.matr = matrix(pratemap,nrow = dims[1], ncol = dims[2]*dims[3])
  pratemap.diff = diff(pratemap.matr)
  dim(pratemap.diff) = dim(pratemap) - c(1,0,0)
  rm(pratemap.matr)
  
  model = list()
  model$cfg = cfg
  
  print("training base RBMs")
  model$rbm.l0 = model.getBaseRbms(model$cfg, cfg$rbm.l0, pratemap.diff,)
  if (andsave)
    model.save(model,model$cfg)

  rbmup.l0 = model.makeRBMup.base(model$rbm.l0,
                                        model$cfg,
                                        model$cfg$rbm.l0,
                                        pratemap.diff)  
  
  print("training serial RBMs")
  model$srbm.l1 = model.getSerialRBMs(model$cfg, model$cfg$srbm.l1, cfg$rbm.l0,rbmup.l0)
  if (andsave)
    model.save(model,model$cfg)
  
  
  print("training ARLMs")
  model$arlms = model.getLinearModels.arlm(cfg,cfg$pred.l0,cfg$rbm.l0,rbmup.l0)
  print("training complete")
  if (andsave)
    model.save(model,model$cfg)
  
  return(model)
}
