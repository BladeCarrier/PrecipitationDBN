#default model, used for comparison
library(foreach)
library(doSNOW)
library(caret)


predict.byPrev = function(precedingMap)
{
  return(precedingMap)
}


prate.transformed = transform(prate.original)
dim(prate.transformed) = dim(prate.original)

#warning: errors are calculated ontransformed values, see distribution
evaluate = function(daysForward)
{#what takes it so long?!
  dataSplit = daysForward+ createDataPartition((1):(dim(prate.transformed)[1]-daysForward), p = 0.2, list = FALSE )
  
  
  errors = vector(mode = "double")
  for (day in dataSplit)
  {
    d = predict.byPrev(prate.transformed[day-daysForward,,])
    errors = append(errors,prate.rmse(d,prate.transformed[day,,]))
  }
  meanError = sum(errors)/length(errors)
  return(meanError)
}
qplot(y = mapply(evaluate,1:30), x = 1:30)

prate.plotMap.fromDay(prate.restore(d))
prate.plotMap.fromDay(prate.original[3,,])
