install.packages(c("caret","ncdf","doSNOW","foreach","deepnet","penalized"),dependencies = TRUE)



source("Rsource/preproc.R")
source("Rsource/model.R")

#todo list:
# check this thing works:
# - symmetric diaps
# - grid RBM maps
# - sequential RBM
#adapt old 100n (and 200n if it's ok) to fit new diaps and grid structure

#0) try l1 rbm over sequence of [order] l0 rbms as prediction source
#1) predict method and getLinearModels method using both state and difference
#2) check model optimal lambda, mb try lower/higher
#3) find best tune for one-layer rbm on CV data: minimize constant zeroes, maximize out-of-sample clusterability (in guide: likelihood)
#4) experiment with time'of'year codes, krigging decay, multilevels, tree predictors, etc.
#make generic in-batch parallel = function(batchsize,taskIndices, innerFunction) to use in makeRBMs, etc?
#mb use lonDiap and latDiap as model cfg parts instead of computing em?

#в трансформе:
#prate.rmse(pe,pb) #эталон-база
#0.02736182
#prate.rmse(pe,pp) #эталон-прогноз (100x16x12 o8x6)
#0.01137628
#0.013 for 200 nodes/layer
#ошибка RBM:
#100 hidden 16/12 8/6off
#prate.rmse(pratemap, RBMdowns) 0.007804673
#вывод - надо повышать точность энкодера и перейти, наконец, к cv
cfg = list()
cfg$numThreads=4
cfg$batchSize = cfg$numThreads
cfg$path = "models/"
cfg$name = "n50xtest +state regL2 0.5"

dims = dim(prate.transformed)


#spatial rbm
cfg$rbm.l0 = list()
cfg$rbm.l0$sector.width = 16
cfg$rbm.l0$sector.height = 12
cfg$rbm.l0$layerLonDiap = model.grid.createAxis(8,dims[2])
cfg$rbm.l0$layerLatDiap = model.grid.createAxis(6,dims[3])
cfg$rbm.l0$numHidden = 50
cfg$rbm.l0$numEpochs = 200
cfg$rbm.l0$learningRate = 0.05
cfg$rbm.l0$learningRateDecay = 0.999
cfg$rbm.l0$encoder = encoder.init(diap = c(-0.75,0.75),numClusters = 15)

#serial rbm on top of spatial l0
cfg$srbm.l1 = list()
cfg$srbm.l1$order = 3
cfg$srbm.l1$sector.width = 3
cfg$srbm.l1$sector.height = 3
cfg$srbm.l1$layerLonDiap = model.grid.createAxis(1,length(cfg$rbm.l0$layerLonDiap))
cfg$srbm.l1$layerLatDiap = model.grid.createAxis(1,length(cfg$rbm.l0$layerLatDiap))
cfg$srbm.l1$numHidden = 50
cfg$srbm.l1$numEpochs = 10
cfg$srbm.l1$learningRate = 0.1
cfg$srbm.l1$learningRateDecay = 0.99

#cfg$rbm.state.l0$encoder = encoder.init(diap = c(-0.15,1.05),numClusters = 12)

#autoregressive model
cfg$pred.l0 = list()
cfg$pred.l0$order = 3
cfg$pred.l0$reg.L2 = 0.5
#kriging decays?


#model benchmark
model = model.assemble(cfg,prate.transformed,andsave = TRUE)
pred = model.predict(model,prate.transformed)
save(pred,file = "data/pred.Rdata")

srbmup = model.makeSRBMup(model$cfg,
                             model$cfg$srbm.l1,
                             cfg$rbm.l0,
                             model$srbm.l1,
                             rbmup.l0)

model.cfg=model$cfg
layer.cfg=model$cfg$srbm.l1
base.cfg=cfg$rbm.l0
srbms=model$srbm.l1
baseRBMup=rbmup.l0
#encoder benchmark
pratemap = prate.transformed
pratemap.matr = matrix(pratemap,nrow = prate.dims[1], ncol = prate.dims[2]*prate.dims[3])
pratemap.diff = diff(pratemap.matr)
dim(pratemap.diff) = dim(pratemap) - c(1,0,0)
rm(pratemap.matr)

rbmup = model.makeRBMup.base(model$rbm.l0,cfg,cfg$rbm.l0,pratemap.diff)
rbmdown = model.makeRBMdown.base(model$rbm.l0,cfg,cfg$rbm.l0,rbmup,dim(pratemap.diff))
save(rbmdown,file = "data/updown.Rdata")

#report:
#enc: 15 clusters (-0.75,0.75). decay 50
#enc/dec rmse: 2.624857*10^-7

#rbm: 50 hidden, 200 epochs,LR 0.1, 10x10 dim, 5x5 offset
#code rmse: 0.009973213 (griat!). Flattens a bit.

#rbm: 50 hidden, LR 0.03, 600 epochs
#rmse(RBMdowns,pratemap) 0.01036832
#prate.rmse(RBMdowns*2,pratemap) 0.009761155

#100 hidden 16/12 8/6off
#prate.rmse(pratemap, RBMdowns) 0.007804673

#100 hidden 10/10 x 5/5
#prate.rmse(RBMdowns,pratemap) 0.006079486

#25 hidden 7x5 o 4x3
#prate.rmse(rbmdown,pratemap.diff) 0.02042611




#chatter:
#when ready with RBMs, only retrain prediction models in experiments, keep rbms when possible
#might want to identify zero-variance images and skip training LMs for them and from them.
#might want a top-down predictions: predicted Li+1 is a parameter for predicted Li

#191 112 639 ~2931996
#srv 447 422 879 ~29031996

#alarma: numEpochs and other stuff for RBMs needs revising when it comes to read prediction tests
#high-level RBMS OVER SEQUENCES OF LWOLV?