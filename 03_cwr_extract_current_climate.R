# CWR pre-breeding characterising testing environments: Extract current climate data taking into account crop cycle
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#    Load AGMERRA data          #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# R options

agmerra_dir <- Sys.info(); agmerra_dir <- agmerra_dir[names(agmerra_dir)=="sysname"]
if(agmerra_dir == "Linux"){
  wk_dir <- '/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/_primary_files'; setwd(wk_dir); rm(wk_dir)
} else {
  if(agmerra_dir == "Windows"){
    wk_dir <- "//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/_primary_files"; setwd(wk_dir); rm(wk_dir)
  }
}; 
rm(agmerra_dir)

current_climate <- function(crop, occ_data){
  
  # Extract TMAX
  tmax <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmax_daily_ts_agmerra_1980_2010.nc")
  tmax <- raster::rotate(tmax)
  tmax_data <- cbind(occ_data, raster::extract(x = tmax, y = occ_data[,c("lon","lat")])); rm(tmax)
  
  # Extract TMIN
  tmin <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmin_daily_ts_agmerra_1980_2010.nc")
  tmin <- raster::rotate(tmin)
  tmin_data <- cbind(occ_data, raster::extract(x = tmin, y = occ_data[,c("lon","lat")])); rm(tmin)
  
  # Extract SRAD
  
  
  # Extract PREC
  
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#             TMAX              #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
##Loop with lapply everybody  years##
agmerra_data_tmax <- list.files(pattern='_tmax*', full.names=TRUE)
agmerra_data_tmax <- agmerra_data_tmax[-1]
year <- 1981:2010
library(parallel)
tabla_tmax <-lapply( 1:2, function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = ""))  
  test_tmax <- raster::brick(agmerra_data_tmax[i], var.name = "tmax")
  test_tmax <- rotate(test_tmax)  
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(test_tmax), 369), vmode="double")
    z <- test_tmax
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:369] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat', as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]] <- temp.m
  }else { # Cases where we have normal years
    temp.m <- ff(1, dim=c(ncell(test_tmax), 368), vmode="double")
    z <- test_tmax
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:368] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat',as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]]<-temp.m
  }
  removeTmpFiles(h=0)
  return(temp.dt)
})
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla_tmax, file=paste(counDir, '/rtmax.rds', sep=''))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#             TMIN              #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

##Loop with lapply everybody  years##
agmerra_data_tmin <- list.files(pattern='_tmin*', full.names=TRUE)
agmerra_data_tmin <- agmerra_data_tmin[-1]
year <- 1981:2010
tabla_tmin <-lapply( 1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = ""))  
  test_tmin <- raster::brick(agmerra_data_tmin[i], var.name = "tmin")
  test_tmin <- rotate(test_tmin)  
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(test_tmin), 369), vmode="double")
    z <- test_tmin
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:369] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat', as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]] <- temp.m
  }else { # Cases where we have normal years
    temp.m <- ff(1, dim=c(ncell(test_tmin), 368), vmode="double")
    z <- test_tmin
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:368] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat',as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]]<-temp.m
  }
  removeTmpFiles(h=0)
  return(temp.dt)
})
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla, file=paste(counDir, '/rtmin.rds', sep=''))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#        RADIACION SOLAR        #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

##Loop with lapply everybody  years##
agmerra_data_srad <- list.files(pattern='_srad*', full.names=TRUE)
agmerra_data_srad <- agmerra_data_srad[-1]
year <- 1981:2010
tabla_srad <-lapply( 1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = ""))  
  test_srad <- raster::brick(agmerra_data_srad[i], var.name = "tmin")
  test_srad <- rotate(test_srad)  
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(test_srad), 369), vmode="double")
    z <- test_srad
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:369] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat', as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]] <- temp.m
  }else { # Cases where we have normal years
    temp.m <- ff(1, dim=c(ncell(test_srad), 368), vmode="double")
    z <- test_srad
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:368] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat',as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]]<-temp.m
  }
  removeTmpFiles(h=0)
  return(temp.dt)
})
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla, file=paste(counDir, '/rsrad.rds', sep=''))

###################################################################################################################################


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
suppressMessages(if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} else {library(doParallel)})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#     Load CHIRPS data          #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)
chirps_dir <- Sys.info(); chirps_dir <- chirps_dir[names(chirps_dir)=="sysname"]
if(chirps_dir == "Linux"){
  wk_dir <- "/mnt/data_cluster_4/observed/gridded_products/chirps/daily"; setwd(wk_dir); rm(wk_dir)
} else {
  if(chirps_dir == "Windows"){
    wk_dir <- "//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily"; setwd(wk_dir); rm(wk_dir)
  }
}; 
rm(chirps_dir)


##Loop with lapply everybody  years##
year <- 1981:2010
chirps_data_prec <- list.files(pattern=".tif$", full.names=TRUE)
system.time(data_year <- mclapply(1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = "")) 
  rasters  <- raster::stack(chirps_data_prec[grep(pattern=year[[i]], x=chirps_data_prec)])
  removeTmpFiles(h=0)
  return(rasters)
}, mc.cores = 20))


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#         Precipitación         #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
##Loop with lapply everybody  years##
data_year
modo <- raster::brick("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmean_daily_ts_agmerra_1980_2010.nc ")
modo <- modo[[1]]
modo <- rotate(modo)
year <- 1981:1982
library(parallel)
tabla_prec <-lapply( 1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = ""))
  test_prec <- raster::stack(data_year[i])
  test_prec <- raster::unstack(test_prec)
  test_prec <- mclapply(X=test_prec, FUN = function(x){raster::resample(x, y = modo, method = "ngb")}, mc.cores = 20)
  system.time(test_prec <-  mclapply( 1:length(year),function(i){
    system.time(test<- raster::resample(test_prec, modo, method = "ngb"))
    removeTmpFiles(h=0)
    return(test)
  },mc.cores = 20))
   test_prec1 <- raster::stack(test_prec)
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(test_prec1), 369), vmode="double")
    z <- test_prec1
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:369] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat', as.character(fecha))
    temp.m <- as.data.frame(temp.m)
   temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]] <- temp.m
  }else { # Cases where we have normal years
    temp.m <- ff(1, dim=c(ncell(test_prec1), 368), vmode="double")
    z <- test_prec1
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:368] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat',as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]]<-temp.m
  }
  removeTmpFiles(h=0)
  return(temp.dt)
}) 
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(test_prec, file=paste(counDir, '/rprec.rds', sep=''))








































##################################################################################################################
##system.time(resample(test_prec,modo ), gcFirst = TRUE)    
##Loop with lapply everybody  years##
data_year
modo <- raster::brick("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmean_daily_ts_agmerra_1980_2010.nc ")
modo <- modo[[1]]
modo <- rotate(modo)
year <- 1981:2010
library(parallel)
tabla_prec <-lapply( 1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = ""))
  test_prec <- raster::stack(data_year[i])
  test_prec <- raster::unstack(test_prec)
  test_prec <- mclapply(test_prec, FUN = function(i){raster::resample(test_prec,  modo, method = "ngb")}, mc.cores = 20)
  test_prec <- raster::stack(test_prec)
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(test_prec), 369), vmode="double")
    z <- test_prec
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:369] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat', as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]] <- temp.m
  }else { # Cases where we have normal years
    temp.m <- ff(1, dim=c(ncell(test_prec), 368), vmode="double")
    z <- test_prec
    t <- getValues(z)
    t[which(t==-9999)] <- NA
    na.omit(t)
    temp.m[,1] <- 1:nrow(t)
    temp.m[,2:3] <- xyFromCell(object=z, cell=1:nrow(t))
    temp.m[,4:368] <- t[]
    temp.m <- as.ffdf(temp.m)
    names(temp.m) <- c('cellID', 'lon', 'lat',as.character(fecha))
    temp.m <- as.data.frame(temp.m)
    temp.m <- data.table(temp.m)
    temp.m <- temp.m[complete.cases(temp.m),]
    temp.dt[[i]]<-temp.m
  }
  removeTmpFiles(h=0)
  return(temp.dt)
}) 
