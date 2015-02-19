model.getBaseRbms = function(model.cfg, layer.cfg, pratemap,andsave = FALSE)
{
  #why not using model as argument: model might be several GB large to flood down the entire process.
  

  numHidden = layer.cfg$numHidden
  numEpochs = layer.cfg$numEpochs
  learnrate = layer.cfg$learningRate
  lrdecay = layer.cfg$learningRateDecay
  sectorWidth = layer.cfg$sector.width
  sectorHeight = layer.cfg$sector.height
  
  modelLonDiap = layer.cfg$layerLonDiap
  modelLatDiap = layer.cfg$layerLatDiap
  getDiap = model.grid.getDiap
  
  rbms = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(rbms) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
  dims = dim(pratemap)
  gc()
  
  cl <- makeCluster(model.cfg$numThreads)
  registerDoSNOW(cl)
  
  encoder = layer.cfg$encoder
  encode = encoder.encode
  for (lon in modelLonDiap)
  {
    print(paste("processing lon",lon))
    lres = foreach(lat=modelLatDiap, .combine='cbind',.noexport = "rbms") %dopar% 
    {
      library(deepnet)
      
      dims = dim(pratemap)
      
      
      lonDiap = getDiap(lon,sectorWidth,dims[2])
      latDiap = getDiap(lat, sectorHeight,dims[3])
      cut = pratemap[,lonDiap,latDiap]
      cut = encode(encoder,cut)
      
      mm = matrix(cut,nrow = dims[[1]],ncol = length(lonDiap)*length(latDiap)*encoder$numClusters)
      model = rbm.train(mm,hidden = numHidden,numepochs = numEpochs,learningrate = learnrate,learningrate_scale = lrdecay)  
      gc()
      c(paste("tile.",lon,'.',lat,sep=""), list(model))
      
    }
    for (r in 1:ncol(lres))
    {
      rbms[[lres[[1,r]]]] = lres[[2,r]]
    }
  }
  stopCluster(cl)
  if(andsave)
    model.save(rbms,model.cfg,mark = ".rbms")
  return(rbms)
}

model.makeRBMup.base = function(rbms,model.cfg, rbm.cfg,pratemap)
{
  require(foreach)
  require(doSNOW)
  require(deepnet)
  
  print("performing parallel RBM up")
  #strategy: split data into batches and execute parallel algos for batches, not to transfer the whole model
  
  
  numHidden = rbm.cfg$rbm.numHidden
  sectorWidth = rbm.cfg$sector.width
  sectorHeight = rbm.cfg$sector.height
  dims =dim(pratemap)  
  
  encoder = rbm.cfg$encoder
  encode = encoder.encode
  
  batchSize = model.cfg$batchSize
  numThreads = model.cfg$numThreads
  
  modelLonDiap = rbm.cfg$layerLonDiap
  modelLatDiap = rbm.cfg$layerLatDiap
  getDiap = model.grid.getDiap
  
  RBMups = array(list(), dim = c(length(modelLonDiap),length(modelLatDiap)))
  names.lookup = expand.grid(modelLatDiap,modelLonDiap)
  names(RBMups) = paste("tile.",names.lookup[,2],".",names.lookup[,1] ,sep = "")
  rm(names.lookup)
  
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
    
    rbmLocal = list()
    for(i in 1:dim(batch)[1]) 
    {
      rbmind = paste("tile.",batch[i,1],'.',batch[i,2],sep="")
      rbmLocal[[rbmind]] = rbms[[rbmind]]
    }
    
    
    lres = foreach(i = 1:dim(batch)[1],.combine = "cbind",.noexport = "model") %dopar%
    {
      library(deepnet)
      lon = batch[i,1]
      lat = batch[i,2]
      rbm = rbmLocal[[paste("tile.",lon,'.',lat,sep="")]]
      lonDiap = getDiap(lon,sectorWidth,dims[2])
      latDiap = getDiap(lat, sectorHeight,dims[3])
      cut = pratemap[,lonDiap,latDiap]
      cut = encode(encoder,cut)
      
      mm = matrix(cut,nrow = dims[[1]],ncol = length(lonDiap)*length(latDiap)*encoder$numClusters)
      up = rbm.up(rbm,mm)
      gc()  
      c(paste("tile.",lon,'.',lat,sep=""), list(up))
      
      
    }
    dim(lres) = c(2,nrow(batch))#eliminating one-element batch case (ncol = NUL
    for (r in 1:ncol(lres))
    {
      RBMups[[lres[[1,r]]]] = lres[[2,r]]
    }
    
    
    if (batchSize >= dim(pairCoords)[1])
      break
    pairCoords = pairCoords[(batchSize+1):dim(pairCoords)[1],]
    
  }
  stopCluster(cl)
    
  gc()
  return(RBMups)

}
model.makeRBMdown.base = function(rbms,model.cfg, rbm.cfg,rbmup,dims)
{
  
  require(foreach)
  require(doSNOW)
  require(deepnet)
  
  values = array(0.,dim = dims)
  normalizers =  array(0.,dim = dims)
  
  print("performing parallel average RBM down")
  #strategy: split data into batches and execute parallel algos for batches, not to transfer the whole model

  numHidden = rbm.cfg$rbm.numHidden
  sectorWidth = rbm.cfg$sector.width
  sectorHeight = rbm.cfg$sector.height
  
  encoder = rbm.cfg$encoder
  encode = encoder.encode
  decode = encoder.decode
  
  batchSize = model.cfg$batchSize
  numThreads = model.cfg$numThreads
  
  modelLonDiap = rbm.cfg$layerLonDiap
  modelLatDiap = rbm.cfg$layerLatDiap
  getDiap = model.grid.getDiap
  
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
    rbmLocal = list()
    for(i in 1:dim(batch)[1]) 
    {
      rbmind = paste("tile.",batch[i,1],'.',batch[i,2],sep="")
      rbmUpLocal[[rbmind]] = rbmup[[rbmind]] 
      rbmLocal[[rbmind]] = rbms[[rbmind]]
    }
    
    
    lres = foreach(i = 1:dim(batch)[1],.combine = "cbind",.noexport = "model") %dopar%
    {
      library(deepnet)
      lon = batch[i,1]
      lat = batch[i,2]
      #check if this works for several days of rbmup at once
      lonDiap = getDiap(lon,sectorWidth,dims[2])
      latDiap = getDiap(lat,sectorHeight,dims[3])
      
      rbm = rbmLocal[[paste("tile.",lon,'.',lat,sep="")]]
      up = rbmUpLocal[[paste("tile.",lon,'.',lat,sep="")]]
      map=rbm.down(rbm,up)
      dim(map) = c(dims[1],length(lonDiap),length(latDiap),encoder$numClusters)
      map=decode(encoder,map)
      dim(map) = c(dims[1],length(lonDiap),length(latDiap))
      #why not using sectorWidth: in age cases, length(lonDiap) < sectorWidth as the sector is on the edge of map
      
      
      gc()  
      c(paste("tile.",lon,'.',lat,sep=""), list(map),list(lonDiap),list(latDiap))
      
      
    }
    dim(lres) = c(4,nrow(batch))#eliminating one-element batch case (ncol = NULL)
    for (r in 1:ncol(lres))
    {
      map = lres[[2,r]]
      lonDiap = lres[[3,r]]
      latDiap = lres[[4,r]]
      
      #chunk used to avoid dimensionality bugs for edge cases where dim[i] = 1
      chunk = values[,lonDiap,latDiap]
      dim(chunk) = dim(map)
      values[,lonDiap,latDiap] = map + chunk
      normalizers[,lonDiap,latDiap] = 1+normalizers[,lonDiap,latDiap]
    }
    
    
    if (batchSize >= dim(pairCoords)[1])
      break
    pairCoords = pairCoords[(batchSize+1):dim(pairCoords)[1],]
    
  }
  stopCluster(cl)
    
  gc()
  ans = (values/normalizers)
  return(ans)
}
