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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#   siembra y cosecha 50 Km   irigado   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
planting_ir_ggcmi <- raster::brick("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_ir_growing_season_dates_v1.25.nc4", varname = "planting day")
planting_ir_ggcmi <- planting_ir_ggcmi[[1]]
harvest_ir_ggcmi <- raster::brick("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_ir_growing_season_dates_v1.25.nc4", varname = "harvest day")
harvest_ir_ggcmi <- harvest_ir_ggcmi[[1]]

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Cambio de resolucion            #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
tmean <- raster::brick("//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/tmean_daily_ts_agmerra_1980_2010.nc ")
tmean <- tmean[[1]]
tmeanr<-rotate(tmean)
planting_ir_ggcmi1 <-resample(planting_ir_ggcmi,tmeanr )
harvest_ir_ggcmi1  <-resample(harvest_ir_ggcmi,tmeanr )

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#  Quitar los valores -99 y negativos   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
planting_ir_ggcmi1[which(planting_ir_ggcmi[] == -99)] <- NA
planting_ir_ggcmi1[which(planting_ir_ggcmi1[] < 0)] <- NA
harvest_ir_ggcmi1[which(harvest_ir_ggcmi[] == -99)] <- NA
harvest_ir_ggcmi1[which(harvest_ir_ggcmi1[] < 0)] <- NA

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#      Ciclo del Cultivo Irrigado       #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
ciclo_ir <- overlay(harvest_ir_ggcmi1,
                 planting_ir_ggcmi1,
                 fun=function(r1,r2){return(abs(r1-r2))})

countries <- rgdal::readOGR(dsn = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/world_shape/world_shape", "all_countries")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
plot(ciclo_ir, main= "Crop Irrigado system")
plot(countries, add = T)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#    siembra y cosecha 50 Km   secano   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
planting_rf_ggcmi <- raster::brick("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
harvest_rf_ggcmi <- raster::brick("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Cambio de resolucion            #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
planting_rf_ggcmi1 <-resample(planting_rf_ggcmi,tmeanr )
harvest_rf_ggcmi1  <-resample(harvest_rf_ggcmi,tmeanr )

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#  Quitar los valores -99 y negativos   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
planting_rf_ggcmi1[which(planting_rf_ggcmi[] == -99)] <- NA
planting_rf_ggcmi1[which(planting_rf_ggcmi1[] < 0)] <- NA
harvest_rf_ggcmi1[which(harvest_rf_ggcmi[] == -99)] <- NA
harvest_rf_ggcmi1[which(harvest_rf_ggcmi1[] < 0)] <- NA

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Ciclo del Cultivo Secano        #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
ciclo_rf <- overlay(harvest_rf_ggcmi1,
                 planting_rf_ggcmi1,
                 fun=function(r1,r2){return(abs(r1-r2))})
countries <- rgdal::readOGR(dsn = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/world_shape/world_shape", "all_countries")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
plot(ciclo_rf, main="Crop Secano system")
plot(countries, add = T)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#    Diferencia irrigado vs secano      #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
ciclo_ir_rf <- overlay(ciclo_ir,
                      ciclo_rf,
                      fun=function(r1,r2){return(r1-r2)})
ciclo_ir_rf[which(ciclo_ir_rf[]<0)]<-NA
countries <- rgdal::readOGR(dsn = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/world_shape/world_shape", "all_countries")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
plot(ciclo_ir_rf,main="Crop Irrigado vs Secano")
plot(countries, add = T)











