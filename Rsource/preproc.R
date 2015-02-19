
library("ncdf")
library("ggplot2")





prate.transform = function(p)
{
  return(0.1+2*log((1+3*p^0.25)))
}
prate.restore = function(p)
{
  return ( ((exp((p-0.1)/2) -1)/3)^4 )
}
prate.rmse = function(map1,map2)
{
  return (sum((map1 - map2)^2)/length(map1))
}

prate.toMap = function( dayMap)
{
  pts = data.frame(1:length(dayMap))
  pts$x=(1:length(dayMap) -1)%%nrow(dayMap)+1
  pts$y=(1:length(dayMap)-1)%/%nrow(dayMap)+1
  pts$prate=matrix(dayMap,nrow =length(dayMap))
  return(pts)
}
prate.plotMap.fromSeries = function(prate,day,title = "plot")
{
  prate.plotMap.fromDay(prate[day,,],title)
}
prate.plotMap.fromDay = function(dayMap,title = "plot")
{
  pts = prate.toMap(dayMap)
  longitude = pts$x * 360/192.
  latitude = pts$y * 180 / 94. - 90
  prate = pts$prate
  qplot( longitude,lattitude, color =  -prate,main = title)  
}
prate.plotMap.fromDay.lattice = function(dayMap,title = "plot")
{
  require(lattice)
  xmin = 0
  xmax = max(dayMap)
  
  levelplot(as.matrix(dayMap),
            main=title,
            xlab="Longitude Units",
            ylab="Latitude Units",
            col.regions=colorRampPalette(c("lightblue", "blue"), space = "rgb")(120),
            cuts=100,
            at=seq(xmin, xmax, (xmax-xmin)/20))
}
nc = open.ncdf("data/98map.nc")
prate = get.var.ncdf(nc,"prate")
prate[prate < 0] = 0
prate = aperm(prate, c(3,1,2))

prate.original = prate
prate.transformed = prate.transform(prate)
prate.dims = dim(prate.transformed)


