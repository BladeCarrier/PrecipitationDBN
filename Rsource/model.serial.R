model.getSerialRBMs = function(model.cfg,layer.cfg,base.cfg,baseRBMup,andsave = FALSE)
{

  
  numThreads = model.cfg$numThreads

  
  order = layer.cfg$order
  numHidden = layer.cfg$numHidden
  numEpochs = layer.cfg$numEpochs
  learnrate = layer.cfg$learningRate
  lrdecay = layer.cfg$learningRateDecay
  sectorWidth = layer.cfg$sector.width
  sectorHeight = layer.cfg$sector.height  
  
  dayDiap = (order):dim(baseRBMup[[1]])[1]
  
  
  modelLonDiap = layer.cfg$layerLonDiap
  modelLatDiap = layer.cfg$layerLatDiap
  
  getDiap = model.grid.getDiap 
  
  srbms = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(srbms) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
  gc()
  dims = dim(baseRBMup)
  
  cl <- makeCluster(numThreads)
  registerDoSNOW(cl)
  print("training LMs")
  for (lon in modelLonDiap)
  {
    print(paste("now processing lon", lon))
    lres = foreach(lat = modelLatDiap,.combine = "cbind") %dopar%
    {
      library(deepnet)
      library(foreach)
      lonDiap = getDiap(lon,sectorWidth,dims[1])
      latDiap = getDiap(lat,sectorHeight,dims[2])
      cut = baseRBMup[lonDiap,latDiap]

      
      mm = as.matrix( foreach(day = dayDiap,.combine = "rbind") %do%
      {#collect all data triples.
        dayFrame = (day - order + 1):day
        fres = foreach(ri = 1:(length(lonDiap)*length(latDiap)),.combine = "cbind") %do%
        {
          cut[[ri]][dayFrame,]
        }#dim = c(dayFrameLength,sum(numHidden))
        as.vector(fres)
      })
      
      model = rbm.train(mm,hidden = numHidden,numepochs = numEpochs,learningrate = learnrate,learningrate_scale = lrdecay)  
      c(paste("tile.",lon,'.',lat,sep=""), list(model))
      
    }
    for (r in 1:ncol(lres))
    {
      srbms[[lres[[1,r]]]] = lres[[2,r]]
    }
  }
  stopCluster(cl)
  gc()
  if(andsave)
    model.save(srbms,model.cfg,mark = ".srbms")
  return(srbms)
  
}



model.makeSRBMup = function(model.cfg,layer.cfg,base.cfg,srbms,baseRBMup,andsave = FALSE)
{
  
  
  numThreads = model.cfg$numThreads
  
  
  order = layer.cfg$order
  numHidden = layer.cfg$numHidden
  numEpochs = layer.cfg$numEpochs
  learnrate = layer.cfg$learningRate
  lrdecay = layer.cfg$learningRateDecay
  sectorWidth = layer.cfg$sector.width
  sectorHeight = layer.cfg$sector.height  
  
  dayDiap = (order):dim(baseRBMup[[1]])[1]
  
  
  modelLonDiap = layer.cfg$layerLonDiap
  modelLatDiap = layer.cfg$layerLatDiap
  
  getDiap = model.grid.getDiap 
  
  RBMups = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(RBMups) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
  gc()
  dims = dim(baseRBMup)
  
  cl <- makeCluster(numThreads)
  registerDoSNOW(cl)
  print("making unitary SRBM up")
  tasks = expand.grid(modelLonDiap,modelLatDiap)
  lres = foreach(i = 1:nrow(tasks),.combine = "cbind") %dopar%
  {
    library(deepnet)
    library(foreach)
    lon =tasks[i,1]
    lat =tasks[i,2]
    lonDiap = getDiap(lon,sectorWidth,dims[1])
    latDiap = getDiap(lat,sectorHeight,dims[2])
    cut = baseRBMup[lonDiap,latDiap]
    srbm = srbms[[paste("tile.",lon,'.',lat,sep="")]]
    
    mm = as.matrix( foreach(day = dayDiap,.combine = "rbind") %do%
    {#collect all data triples.
      dayFrame = (day - order + 1):day
      fres = foreach(ri = 1:(length(lonDiap)*length(latDiap)),.combine = "cbind") %do%
      {
        cut[[ri]][dayFrame,]
      }#dim = c(dayFrameLength,sum(numHidden))
      as.vector(fres)
    })
  
    srbm.up = rbm.up(srbm,mm)
    gc()
    c(paste("tile.",lon,'.',lat,sep=""), list(srbm.up))
  
  }
  for (r in 1:ncol(lres))
  {
    RBMups[[lres[[1,r]]]] = matrix(lres[[2,r]],nrow =length(dayDiap),ncol =numHidden)
  }
  
  stopCluster(cl)
  gc()
  if(andsave)
    model.save(srbms,model.cfg,mark = ".srbms")
  return(srbms)

}