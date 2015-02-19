rbm.toXData = function(numHidden,up,order,day)
{
  entry = vector()
  for (w in 1:numHidden)
  {
    for (t in 1:order)
      entry[paste("x",-t,w)] = up[day-t,w]
  }
  entry
}
rbm.toYData = function(numHidden,up,day)
{
  entry = vector()
  for (w in 1:numHidden)
  {
    entry[[paste("y",w)]] = up[day,w]
  }
  entry
}

linearModel.predict= function(modelMatr, entry)
{
  require(foreach)
  ans = foreach (yimg = 1:dim(modelMatr)[1], .combine = "c") %do%
  {
    sum(modelMatr[yimg,] * c("(Intercept)"=1,entry))
  }
  ans = linearModel.sigma(ans)
  return(ans)
}
linearModel.sigma = function(value)
{
  1/(1 + exp(-value))
}
linearModel.logit = function(value)
{#inverse sigma
  log(value/(1-value))
}

encoder.init = function(numClusters = 15, diap = c(0,1))
{
  #warning: pls have some spare diap from left and right
  encoder = list()
  encoder$numClusters=numClusters
  encoder$diap = diap
  encoder$centers = ((1:encoder$numClusters) - 0.5)/encoder$numClusters
  encoder$centers = encoder$centers * (diap[2]-diap[1]) + diap[1]
  return(encoder)
}

encoder.encode = function(encoder,value,gammaDecay = 50,cutoff = 5)
{
  diap = encoder$diap
  index = round(encoder$numClusters* (value-diap[1])/(diap[2]-diap[1]))
  index = ifelse(index >encoder$numClusters, encoder$numClusters,
                 ifelse(index < 1, 1, index))
  code = array(0, dim = c(length(value),encoder$numClusters))
  for (vi in 1:length(value))
  {
    active = max(index[vi]-cutoff,1):min(index[vi]+cutoff, encoder$numClusters)
    distance = abs(value[vi] - encoder$centers[active])
    code[vi,active] = exp(-gammaDecay*(diap[2]-diap[1])*distance^2)
    #gaussian decay
    
    #normalize; 
    #most recent change; remove comment if it works by 05.02.14
    
  }
  if(!is.null(dim(value)))
  {
    dim(code) = c(dim(value),encoder$numClusters)
  }else
  {
    dim(code) = c(length(value ), encoder$numClusters)
  }
  return (code)
  
}

model.grid.createAxis = function(offset,length,centralize = TRUE)
{
  diap = 0:((length-1)/offset)*offset+1
  if(centralize)
  {
    dist = length - diap[length(diap)]
    diap = diap + dist%/%2
  }
  diap
}
model.grid.getDiap = function(coord, breadth, limit)
{#cover area of one individual rbm/srbm
  leftBound = (coord - (breadth%/%2))
  diap = leftBound:(leftBound + breadth - 1)
  diap = diap[diap > 0 & diap <= limit]
  diap
}

encoder.decode = function(encoder,code,threshold = 0.)
{
  require(foreach)
  dim(code) = c(length(code)/encoder$numClusters, encoder$numClusters )
  code = code*(code>threshold)
  ans = foreach(i = 1:dim(code)[1], .combine = "rbind" ) %do%
  {
    ifelse(sum(code[i,] != 0), 
           weighted.mean(encoder$centers,w= code[i,]),
           encoder$centers[which.max(code)])#if nothing is above threshold, return most likely one
    
    
           
  }
  dim(ans) = dim(code)[1:(length(dim(code))-1)]
  ans
}
