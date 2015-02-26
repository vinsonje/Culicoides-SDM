#########################
#Code for running BRT with Cuilc. Data
#########################
#ADDED NEW CODE TO TEST GITHUB
#Set working directory
setwd("~/Dropbox/Culicoides/Data")
#Save default graph par
graphics.off()
#Load required packages
library(dismo)
library(gbm)
library(ROCR)

set.seed(101)

#Read in data and remove a subset
data.full=read.csv("PA_data.csv")

set.aside = sample(1:dim(data.full)[1], dim(data.full)[1]*0.4, replace=F)
data.set.aside = data.full[set.aside,]

data = data.full[-set.aside,]


#Deciding which species to focus on
#Using GLM
test = glm(data[,77]~data[,2]+data[,3]+data[,4]+data[,5]+data[,6]+data[,7]+data[,8]+data[,9]+data[,10]+data[,11]
           +data[,12]+data[,13]+data[,14]+data[,15]+data[,16]+data[,17]+data[,18]+data[,19]+data[,20]
           +data[,21]+data[,22]+data[,23]+data[,24]+data[,25]+data[,26]+data[,27]+data[,28]+data[,29]+data[,30]
           +data[,31]+data[,32]+data[,33]+data[,34]+data[,35]+data[,36]+data[,37]+data[,38]+data[,39]+data[,40]
           +data[,41]+data[,43]+data[,44]+data[,45]+data[,46]+data[,47]+data[,48]+data[,49]+data[,50]
           +data[,51]+data[,52]+data[,53],
           family=poisson)

#Using the results create boxplots of the disease score 
#only selected those that had a positive coef.
spec.vec = c(20,22,24,3,47)
par(mfrow=c(2,3))
for(i in spec.vec){
  id0 = which(data[,i]==0)
  id1 = which(data[,i]==1)
  
  boxplot(data[id0,77],data[id1,77], main=names(data)[i], names=c("Abs", "Pres"), ylab="Disease Score")
}

par(mfrow=c(2,3))
#Use BRT to generate 5 different models for the species decided from the GLM
#Save these models in a list 
BRT.data = vector(mode="list", length=length(spec.vec))
for(i in 1:length(spec.vec)){
  BRT.data[[i]] = gbm.step(data=data, 
                          gbm.x = 56:74,
                          gbm.y = spec.vec[i],
                          family = "bernoulli",
                          tree.complexity = 5,
                          learning.rate = 0.001,
                          bag.fraction = 0.5)
}

#Order alphabetically so all rel. inf. values are on the same
#row.
BRT.data.sum.ord = vector(mode="list", length=length(spec.vec))
for(i in 1:length(spec.vec)){
  BRT.data_temp = as.data.frame(summary(BRT.data[[i]], plotit=F))
  BRT.data.sum.ord[[i]] = BRT.data_temp[order(BRT.data_temp$var),]
}
#Add up all rel. inf. for each predictor
rel.inf.sum = rep(0, length(56:74))
for(i in 1:length(rel.inf.sum)){
  row.sum=0
    for(j in 1:length(spec.vec)){
      row.sum = row.sum + BRT.data.sum.ord[[j]]$rel.inf[i]
    }
  rel.inf.sum[i]=row.sum
}
#Create dataframe for the added up rel. inf.
data.rel.inf.sum = data.frame(cbind(levels(BRT.data.sum.ord[[1]]$var),rel.inf.sum))
#Figure out the order from smallest to largest inf.
order = order(data.rel.inf.sum$rel.inf.sum)
order = order[1:7]
#Create barplots of the relative influence of each predictor
#Plot in order of TOTAL rel. inf.
par(mar=c(5,12,4,2),las=2, mfrow=c(2,3))
for(j in 1:length(spec.vec)){
BRT.data_test = as.data.frame(summary(BRT.data[[j]], plotit=F))
BRT.data_test = BRT.data_test[order(BRT.data_test$var),]
barplot(BRT.data_test$rel.inf[order], names.arg=BRT.data_test$var[order], horiz=T, main=names(data)[spec.vec[j]])
}


#Plotting ROC curves
par(mfrow=c(2,3), mar=c(5.1, 4.1, 4.1, 2.1))
auc.values = rep(0, length(spec.vec))
for(i in 1:length(spec.vec)){
  test.pred = predict.gbm(BRT.data[[i]], data.set.aside, n.trees=BRT.data[[i]]$gbm.call$best.trees, type='response')
  pred = prediction(test.pred, data.set.aside[,spec.vec[i]])
  perf = performance(pred, "tpr", "fpr")
  auc = performance(pred, "auc")
  auc = as.numeric(auc@y.values)
  auc.values[i] = auc
  plot(perf, xaxs='i', yaxs='i', main=names(data)[spec.vec[i]])
  abline(0,1.0)
}

#####################################
#Species Predictions
#####################################

#load in required packages
#Also need to change the working directory
require(maptools)
require(maps)
require(mapdata)
require(raster)
setwd("~/Dropbox/Culicoides")
par(mfrow=c(1,1))

#Pulling in the data
#Then change the ranges that we are interested in investigating
#Then make it a raster brick.
BClim=getData('worldclim', var='bio', res=2.5, path='data/')
YbrevRange = extent(-95,-70,24.5,37)
BClim = crop(BClim, YbrevRange)
writeRaster(BClim, filename="data/YbrevBC_2.5.grd", overwrite=T)
BClim = brick("data/YbrevBC_2.5.grd")
names(BClim) = names(data)[56:74]

#This plots one environmental variable and the study sites
# plot(BClim, 16, cex=0.5, legend=T, mar=par("mar"), xaxt="n", yaxt="n", main=names(BClim)[16])
# map("state", xlim=c(-95, -70), ylim=c(24.5, 37), fill=F, col="cornsilk", add=T)
# points(data$lon, data$lat, pch=20, cex=0.5, col="darkgreen")
# axis(1,las=1)
# axis(2,las=1)
# box()

#Use the raster brick to make species range predictions
#Also plot the sites that had the species present.
pred.vec = vector(mode="list", length=length(spec.vec))
par(mfrow=c(3,2), mar=c(5, 4, 4, 2))
for(i in 1:length(spec.vec)){
  p = predict(BClim, BRT.data[[i]],n.trees=BRT.data[[i]]$gbm.call$best.trees, type="response")
  pred.vec[[i]]=p

  plot(p, main=names(data)[spec.vec[i]], xlim=c(-95,-75), ylim=c(24.5,37))
  map("state", xlim=c(-95, -75), ylim=c(24.5, 37), fill=F, col="black", add=T)
}

# #####################################
# #Average Disease Score Map
# #####################################
# par(mfrow=c(1,1))
# #Calculate the average disease score of where the species were present
# #of each species we're interested in.
# dis.scr.avgs = rep(0,length(spec.vec))
# for(i in 1:length(spec.vec))
# {
#   pres = which(data[,spec.vec[i]]==1)
#   dis.scr.avgs[i] = mean(data$heat[pres])
# }
# 
# #Create a matrix to have the predicted average disease score
# #of the US
# coor = coordinates(pred.vec[[i]])
# dis.score.map.data = matrix(0, length(coor[,1]), length(spec.vec)+3)
# dis.score.map.data[,1:2] = coor
# 
# #Fill in the matrix with average disease score
# #of the average DS of a species
# #You can change the threshold prediction probability
# #which will change the risk map
# for(i in 1:length(spec.vec))
# {
#   p=pred.vec[[i]]
#   pred.pres = which(values(p)>0.6)
#   dis.score.map.data[pred.pres,i+2] = dis.scr.avgs[i]
# }
# 
# #Create a raster file that has the disease score
# for(i in 1:length(coor[,1])){
#   dis.score.map.data[i,8] = mean(dis.score.map.data[i,3:7])
# }
# pred.dis.score.map = dis.score.map.data[,c(-3,-4,-5,-6,-7)]
# colnames(pred.dis.score.map) = c("x", 'y', 'heat')
# pred.dis.score.map = rasterFromXYZ(pred.dis.score.map)
# 
# #Plot the disease map
# plot(pred.dis.score.map, main="Disease Risk")
# map("state", xlim=c(-95, -75), ylim=c(24.5, 37), fill=F, col="black", add=T)
# 
# #Add points of the sites with the highest risk
# for(i in 1:length(set.aside)){
# lat = data.full$lat[set.aside[i]]
# lon = data.full$lon[set.aside[i]]
# score = as.character(data.full$heat[set.aside[i]])
# points(lon, lat, col="black", pch=score)}
# 
# 
# #############
# #Species Number Map
# #############
# #Create a matrix to have the predicted species number
# #of the US
# spc.score.map.data = matrix(0, length(coor[,1]), length(spec.vec)+3)
# spc.score.map.data[,1:2] = coor
# 
# #Fill in the matrix with species number
# #You can change the threshold prediction probability
# #which will change the risk map
# for(i in 1:length(spec.vec))
# {
#   p=pred.vec[[i]]
#   pred.pres = which(values(p)>0.6)
#   spc.score.map.data[pred.pres,i+2] = 1
# }
# 
# #Create a raster file that has the species number
# for(i in 1:length(coor[,1])){
#   spc.score.map.data[i,8] = sum(spc.score.map.data[i,3:7])
# }
# pred.spc.score.map = spc.score.map.data[,c(-3,-4,-5,-6,-7)]
# colnames(pred.spc.score.map) = c("x", 'y', 'spc score')
# pred.spc.score.map = rasterFromXYZ(pred.spc.score.map)
# 
# #Plot the disease map
# plot(pred.spc.score.map, main="Species Number")
# map("state", xlim=c(-95, -75), ylim=c(24.5, 37), fill=F, col="black", add=T)
# 
# 
# 
# test = gbm.step(data=data, 
#                          gbm.x = 2:53,
#                          gbm.y = 77,
#                          family = "poisson",
#                          tree.complexity = 3,
#                          learning.rate = 0.01,
#                          bag.fraction = 0.5)