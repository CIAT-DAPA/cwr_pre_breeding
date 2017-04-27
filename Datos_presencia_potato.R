# CWR pre-breeding characterising testing environments: base code
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
##   Extraction of raster information   ##
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Load packages           #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
options(warn = -1); options(scipen = 999); g <- gc(reset = T); rm(list = ls())
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf)){install.packages('ncdf'); library(ncdf)} else {library(ncdf)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(doMC)){install.packages('doMC'); library(doMC)} else {library(doMC)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#   Load Monfread and Mapspam   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

mapspam<-raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/Datos_presencia/potato_mapspam.nc")
monfread <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/Datos_presencia/potato_monfread.nc")
na.omit(mapspam[])
na.omit(monfread[])

#### Datos de presencia mapspam
mapspam_threshold <- mapspam
mapspam_threshold[which(mapspam_threshold[] <= 0)] <- NA
mapspam_threshold[which(mapspam_threshold[] > 0)] <- 1

#### Datos de presencia monfread
monfread_threshold <- monfread
monfread_threshold[which(monfread_threshold[] <= 0)] <- NA
monfread_threshold[which(monfread_threshold[] > 0)] <- 1

### Presencia Monfread y Mapspam 
rasterStack <- raster::stack(mapspam_threshold, monfread_threshold)
rasterSum <- sum(rasterStack, na.rm = T)

### Datos de Genesys  Common crop
datos <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/Datos_presencia/genesys/geo.csv", header = T)
l1<-datos[,2]
l2<-datos[,3]
### Datos de Genesys Wild and others
dat <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/Datos_presencia/genesys/sil.csv", header = T)
l3<-dat[,2]
l4<-dat[,3]
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#             Resultados               #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-#

countries <- rgdal::readOGR(dsn = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/world_shape/world_shape", "all_countries")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
plot(rasterSum, main="Datos de Presencia Monfread, MapSpam, Genesys")
plot(countries, add = T)
(points(x=l2,y=l1,col=2,pch=20)) # Common crop
(points(x=l4,y=l3,col=4,pch=21)) # Wild and others


