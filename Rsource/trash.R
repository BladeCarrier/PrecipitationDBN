model.predictNext = function(model,lastmap)
{
  #lastmap: several steps of prate map, last one being most recent.
  predicted = array(0.,dim = dim(lastmap[1,,]))
  normalizers =  array(0.,dim = dim(lastmap[1,,]))
  
  rbms = model$rbms
  arlms = model$arlms
  cfg = model$cfg
  
  offset.h = cfg$sector.offset.horizontal.l0
  offset.v = cfg$sector.offset.vertical.l0
  numHidden = cfg$rbm.numHidden.l0
  sectorWidth = cfg$sector.width.l0
  sectorHeight = cfg$sector.height.l0
  order = cfg$model.order.l0
  wdir = cfg$model.path
  dims =dim(lastmap)  
  
  toXData = rbm.toXData
  encoder = model$cfg$encoder.l0
  
  modelLonDiap = 0:((dims[2]-1)/offset.h )*offset.h+1
  modelLatDiap = 0:((dims[3]-1)/offset.v)*offset.v+1
  for(lon in modelLonDiap)
    for(lat in modelLatDiap)
    {
      rbm = rbms[[paste("rbm.",lon,'.',lat,sep="")]]
      arlm = arlms[[paste("arlm.",lon,'.',lat,sep="")]]
      
      lonDiap = lon:min(dims[2],lon+sectorWidth-1)
      latDiap = lat:min(dims[3],lat+sectorHeight-1)
      cut = lastmap[,lonDiap,latDiap]
      cut =  encoder.encode(encoder,cut)
      mm = as.matrix(x = cut)
      dim(mm) = c(dim(cut)[1],length(lonDiap)*length(latDiap)*encoder$numClusters)
      up = rbm.up(rbm,mm)
      
      datarow = toXData(rbm$size,up,order,day=dims[1]+1)
      output = linearModel.predict(arlm,datarow)
      
      map=rbm.down(rbm,output)
      map=encoder.decode(encoder,map)#try apply me
      dim(map) = c(length(lonDiap),length(latDiap))
      predicted[lonDiap,latDiap] = map+predicted[lonDiap,latDiap]
      normalizers[lonDiap,latDiap] = 1+normalizers[lonDiap,latDiap]
      
    }
  ans = (predicted/normalizers)
  return(ans)
}
model.predictAll = function(model, pratemap)
{
  require(foreach)
  require(doSNOW)
  cfg = model$cfg
  #cl <- makeCluster(cfg$numThreads)
  #registerDoSNOW(cl)
  predNext = model.predictNext
  dayDiap=(cfg$model.order.l0+1):dim(pratemap)[1]
  order = cfg$model.order.l0
  stop("first implement parallel prediction to get rbmup, than do parallel rbm down")
  
  lres = foreach(day = dayDiap,.combine="c") %do%  #would need to transfer whole model
{
  source('Rsource/auxilary.R')
  require(deepnet)
  predNext(model,pratemap[(day - order):(day-1),,])
  
}
#stopCluster(cl)

dims = dim(pratemap)
dim(lres) = c(dims[2],dims[3],dims[1] - order)
lres = aperm(lres,c(3,1,2))
result = pratemap*NA
result[dayDiap,,]=lres
return (result)

}
#obsolete
nnmodel.predict= function(models,rbms,prev.last3)
{
  predicted = array(0.,dim = dim(prev.last3[1,,]))
  normalizers =  array(0.,dim = dim(prev.last3[1,,]))
  
  prev.transformed = transform(prev.last3)
  order = cfg$odel.order.l0
  
  dims =dim(prev.last3)
  modelLonDiap = 0:((dims[2]-1)/3 )*3+1
  modelLatDiap = 0:((dims[3]-1)/3 )*3+1
  for(lon in modelLonDiap)
    for(lat in modelLatDiap)
    {
      lonDiap = lon:min(dims[2],lon+10-1)
      latDiap = lat:min(dims[3],lat+10-1)
      rbm = rbms[[paste("rbm.",lon,'.',lat,sep="")]]
      model = models[[paste("arlm.",lon,'.',lat,sep = "")]]
      cut = prev.last3[,lonDiap,latDiap]
      
      mm = as.matrix(x = cut)
      dim(mm) = c(dim(cut)[1],length(lonDiap)*length(latDiap))
      up = rbm.up(rbm,mm)
      
      entry = vector()
      for (w in 1:rbm$size[2])
      {
        for (t in 1:order)
          entry[[paste("x",-t,w)]] = up[3-t+1,w]
        
      }
      mm = as.matrix(entry)
      dim(mm) = c(1,20*order)
      mm[1,] = entry
      output = nn.predict(model,mm)[1,]
      
      map=rbm.down(rbm,output)
      dim(map) = c(length(lonDiap),length(latDiap))
      predicted[lonDiap,latDiap] = map+predicted[lonDiap,latDiap]
      normalizers[lonDiap,latDiap] = 1+normalizers[lonDiap,latDiap]
    }
  ans = prate.restore(predicted/normalizers)
  return(ans)
  
}


#####
model.getAllRbms = function(model, pratemap,andsave = FALSE)
{
  
  
  cfg = model$cfg
  offset.h = cfg$sector.offset.horizontal.l0
  offset.v = cfg$sector.offset.vertical.l0
  numHidden = cfg$rbm.numHidden.l0
  numEpochs = cfg$rbm.numEpochs.l0
  sectorWidth = cfg$sector.width.l0
  sectorHeight = cfg$sector.height.l0
  
  rbms = list()
  dims = dim(pratemap)
  
  cl <- makeCluster(cfg$numThreads)
  registerDoSNOW(cl)
  
  getRBM = rbm.get
  for (lon in 0:((dims[2]-1)/offset.h )*offset.h+1)
  {
    lres = foreach(lat=0:((dims[3]-1)/offset.v)*offset.v+1, .combine='cbind') %dopar% 
{
  library(deepnet)
  
  model = getRBM(pratemap,lon,lat,sectorWidth,sectorHeight,numHidden,numEpochs)
  c(paste("rbm.",lon,'.',lat,sep=""), list(model))
  
}
for (r in 1:ncol(lres))
{
  rbms[lres[[1,r]]] = list(lres[[2,r]])
}
  }

stopCluster(cl)
if (andsave)
  save(rbms,file = paste(cfg$model.path,"rbms",".RData",sep = ""))
return(rbms)
}

#mark
#save and load/assemble
saveLinearodels = function(pratemap,rbms)
{
  
  offset.h = cfg.sector.offset.horizontal.l0
  offset.v = cfg.sector.offset.vertical.l0
  numHidden = cfg.rbm.numHidden.l0
  sectorWidth = cfg.sector.width.l0
  sectorHeight = cfg.sector.height.l0
  order = cfg.model.order.l0
  wdir = cfg.model.path
  dims =prate.dims  
  modelLonDiap = 0:((dims[2]-1)/offset.h )*offset.h+1
  modelLatDiap = 0:((dims[3]-1)/offset.v)*offset.v+1
  
  cl <- makeCluster(cfg.numThreads)
  registerDoSNOW(cl)
  
  for (lon in modelLonDiap)
  {
    print(paste("now processing lon", lon))
    lres = foreach(lat = modelLatDiap,.combine = "c") %dopar%
{
  library(deepnet)
  library(caret)
  library(foreach)
  rbm = rbms[[paste("rbm",lon,lat)]]
  lonDiap = lon:min(dims[2],lon+sectorWidth-1)
  latDiap = lat:min(dims[3],lat+sectorHeight-1)
  cut = pratemap[,lonDiap,latDiap]
  mm = as.matrix(x = cut)
  dim(mm) = c(dim(cut)[1],length(lonDiap)*length(latDiap))
  up = rbm.up(rbm,mm)
  X=(order+1):dim(up)[1]
  xdata = foreach(day = X, .combine = "rbind") %do%
{
  entry = vector()
  for (w in 1:rbm$size[2])
  {
    for (t in 1:order)
      entry[[paste("x",-t,w)]] = up[day-t,w]
    
  }
  entry
}
xdata = as.data.frame(xdata) 

ydata = foreach(day = X, .combine = "rbind") %do%
{
  entry = vector()
  for (w in 1:rbm$size[2])
  {
    entry[[paste("y",w)]] = up[day,w]
  }
  entry
}
ydata = as.data.frame(ydata)

modeldict = list()
for (yimg in 1:ncol(ydata))
{
  model = train(y = ydata[,paste("y",yimg)],x = xdata,method = "lm")
  modeldict[[paste("model",lon,lat,"img",yimg)]] = model
}
save(modeldict,file = paste(wdir,"model ",lon,",",lat,".RData",sep = ""))
remove(modeldict)
remove(xdata)
remove(ydata)
remove(cut)
remove(mm)
remove(rbm)
c("ans")
}
remove(lres)
  }
stopCluster(cl)

}
assembleLinearModels = function(andsave = FALSE)
{#debug
  #load saved lm files and extract only the weights
  #no need to write it in parallel for the main bottleneck here is disk access time
  offset.h = cfg.sector.offset.horizontal.l0
  offset.v = cfg.sector.offset.vertical.l0
  numHidden = cfg.rbm.numHidden.l0
  sectorWidth = cfg.sector.width.l0
  sectorHeight = cfg.sector.height.l0
  order = cfg.model.order.l0
  wdir = cfg.model.path
  dims = prate.dims
  
  modelLonDiap = 0:((dims[2]-1)/offset.h )*offset.h+1
  modelLatDiap = 0:((dims[3]-1)/offset.v)*offset.v+1
  
  linearModelMatrices = list()
  for (lon in modelLonDiap)
    for (lat in modelLatDiap)
    {
      lonDiap = lon:min(dims[2],lon+sectorWidth-1)
      latDiap = lat:min(dims[3],lat+sectorHeight-1)
      #load modeldict
      load(file = paste(wdir,"model ",lon,",",lat,".RData",sep = ""))
      cMatr = foreach(yimg = 1:numHidden, .combine = "rbind") %do%
{
  model = modeldict[[paste("model",lon,lat,"img",yimg)]]
  cVec = model$finalModel$coefficients
  rm(model)
  cVec
}
linearModelMatrices[[paste("model",lon,lat)]] = as.matrix(cMatr)
rm(modeldict)
gc()
    }

if (andsave)
  save(linearModelMatrices,file = paste(cfg.model.path,"lmMatrices",".RData",sep = ""))

return(linearModelMatrices)
}


#

library(caret)
library(foreach)


####
sampleSize = c(5,5)
library(deepnet)
fullmap = matrix(nrow=90*190*365,ncol=sampleSize[1]*sampleSize[2])

itr = 1
for(lo in 1:190)
  for(la in 1:90)
    for (day in 1:365)
    {
      fullmap[itr,] = as.vector(prate.transformed[day,1:sampleSize[1],1:sampleSize[2]])
      itr = itr+1
    }


model = rbm.train(fullmap,hidden = 20,numepochs = 10)  


#old lr model
dims = dim(prate_adj)
for (lon in 1:dims[2])
  for (lat in 1:dims[3])
  {
    data = preproc(prate_adj,lon,lat)
    trs = data[dataSplit,]#TRainingSet
    tss = data[-dataSplit,]#TeStSet
    sLR = train(targ~., data = trs, method = "lm",trControl = cv)
    tss$targ = predict(sLR, newdata = tss)
    for (ts in 1:length(tss))
    {
      ref[tss$day,lon,lat] = tss$targ
    }
    #sLR = train(targ~., data = trs, method = "bstLs",trControl = cv)
    #sLR = train(targ~., data = trs, method = "rpart",trControl = cv)
    #sLR = train(targ~., data = trs, method = "Boruta",trControl = cv)
    
  }
preproc = function(pratemap,lon,lat)
{
  dims = dim(pratemap)
  
  sectorWidth = 5 #11x11 sector
  sectorHeight = 5 #lon and lat diap halves ~ margins
  
  LonDiap = max(lon-sectorWidth,1):min(lon+sectorWidth-1,dims[2])
  LatDiap = max(lat-sectorHeight,1):min(lat+sectorHeight-1,dims[3])
  
  pdata = data.frame(list(day=1:LEN,targ = 0*(1:LEN)))  
  for (i in LonDiap)
  {
    for (j in LatDiap)
    {
      pdata[paste("prate -1",i,j)] = pratemap[1:LEN,i,j]
    }
  }
  
  
  for (t in 1:(LEN-1))
  {    
    pdata[t,]$targ = as.double(pdata[t+1,][paste("prate -1",lon,lat)])
  }
  
  pdata = pdata[1:LEN-1,]
  
  return (pdata)
}
#t/t;cv
dataSplit = createDataPartition(1:364, p = 0.8, list = FALSE )
cv <- trainControl(method = "folds", p = 0.85)

#dnn plain



library(deepnet)
nmap = prate.transformed
xm = matrix(nmap[1:LEN-1,,],nrow = LEN-1)
ym = matrix(nmap[2:LEN,,],nrow = LEN-1) 
dnn <- dbn.dnn.train(xm,ym,hidden = c(1500,1000,500))

a = prate.restore(as.array(nn.predict(dnn,xm)))
dim(a) =  c(364,192,94)
