# CWR pre-breeding characterising testing environments: Identify crop cycle per pixel
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Load packages           #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
# suppressMessages(if(!require(doMC)){install.packages('doMC'); library(doMC)} else {library(doMC)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(stringr)){install.packages('stringr'); library(stringr)} else {library(stringr)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(corrplot)){install.packages('corrplot'); library(corrplot)} else {library(corrplot)})
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(factoextra)){install.packages('factoextra'); library(factoextra)} else {library(factoextra)})
suppressMessages(if(!require(leaflet)){install.packages('leaflet'); library(leaflet)} else {library(leaflet)})



# Load data
occ_data2 <- read.csv("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/database/potato_genesys.csv")
r <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/database/world_body_waters_2-5.asc")
crs(r)<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Validacion
nrow(occ_data2) ## numero de coordenadas 4412
ex <- raster::extract(x=r,y=occ_data2[,c("longitude","latitude")])
cv <- na.omit(ex)
sum(cv)## numero de coordenadas sum: 4178

error <- (1 -sum(cv)/nrow(occ_data2))*100


coor<- cbind(ex, occ_data2[,c('longitude', 'latitude')])
coorval<- coor[complete.cases(coor),]

## Plots

par(mfrow = c(1,2))
plot(r, main= 'Data from Genesys before validation')
points(x=occ_data2[,"longitude"],y=occ_data2[,"latitude"],col=2,pch=20)

plot (r, main ='Data from Genesys after  validation')
points(x=coorval[,"longitude"],y=coorval[,"latitude"],col=2,pch=20)
