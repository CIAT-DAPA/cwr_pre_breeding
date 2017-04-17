#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
## Extraccion de informacion en tablas  ##
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Load packages           #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
options(warn=-1); options(scipen = 999)
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
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#    Load AGMERRA data          #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)
agmerra_dir <- Sys.info(); agmerra_dir <- agmerra_dir[names(amgerra_dir)=="sysname"]
if(agmerra_dir == "Linux"){
  wk_dir <- '/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/_primary_files'; setwd(wk_dir); rm(wk_dir)
} else {
  if(agmerra_dir == "Windows"){
    wk_dir <- "//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/_primary_files"; setwd(wk_dir); rm(wk_dir)
  }
}; 
rm(agmerra_dir)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#             TMAX              #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
##Loop with lapply everybody  years##
agmerra_data_tmax <- list.files(pattern='_tmax*', full.names=TRUE)
lista_tmax <- lapply( 2:length(agmerra_data_tmax), function(i){ 
  cat(paste("Procesando year:", i, "\n", sep = ""))
  test <- raster::brick(agmerra_data_tmax[i], var.name = "tmax")
  test <- rotate(test)
  return(test)
  removeTmpFiles(h=0)
})


## construccion de tablas t_max
year <- seq(from=1980, to= 2010)
tabla_tmax <- lapply(1:2, function(i){ # 1:length(lista_tmax)
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(lista_tmax[[i]]), 369), vmode="double")
    z <- lista_tmax[[i]]
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
    temp.m <- ff(1, dim=c(ncell(lista_tmax[[i]]), 368), vmode="double")
    z <- lista_tmax[[i]]
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
  return(temp.dt)
  removeTmpFiles(h=0)
})

counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla, file=paste(counDir, '/rtmax.rds', sep=''))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#             TMIN              #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
##Loop with lapply everybody  years##
agmerra_data_tmin <- list.files(pattern='_tmin*', full.names=TRUE)
lista_tmin <- lapply(2:length(agmerra_data_tmin), function(i){ 
  cat(paste("Procesando year:", i, "\n", sep = ""))
  test <- raster::brick(agmerra_data_tmin[i], var.name = "tmin")
  test <- rotate(test)
  return(test)
  removeTmpFiles(h=0)
})


## construccion de tablas t_min
year <- seq(from=1980, to= 2010)
tabla_tmin <- lapply(1:2, function(i){ # 1:length(lista_tmin)
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(lista_tmin[[i]]), 369), vmode="double")
    z <- lista_tmin[[i]]
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
    temp.m <- ff(1, dim=c(ncell(lista_tmin[[i]]), 368), vmode="double")
    z <- lista_tmin[[i]]
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
  return(temp.dt)
  removeTmpFiles(h=0)
})

counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla, file=paste(counDir, '/rtmin.rds', sep=''))




#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#        RADIACION sOLAR        #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
##Loop with lapply everybody  years##
agmerra_data_srad <- list.files(pattern='_srad*', full.names=TRUE)
lista_srad <- lapply(2:length(agmerra_data_srad), function(i){ 
  cat(paste("Procesando year:", i, "\n", sep = ""))
  test <- raster::brick(agmerra_data_srad[i], var.name = "srad")
  test <- rotate(test)
  return(test)
  removeTmpFiles(h=0)
})


## construccion de tablas s_rad
year <- seq(from=1980, to= 2010)
tabla_srad <- lapply(1:2, function(i){ # 1:length(lista_srad)
  is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }
  LeapYear <- is.leapyear(year[i]) 
  temp.dt <- list()
  fecha<-seq(as.Date(paste(year[i], "/01/01", sep = "")), as.Date(paste(year[i], "/12/31", sep = "")),by="day")
  if (LeapYear==TRUE){
    temp.m <- ff(1, dim=c(ncell(lista_srad[[i]]), 369), vmode="double")
    z <- lista_srad[[i]]
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
    temp.m <- ff(1, dim=c(ncell(lista_srad[[i]]), 368), vmode="double")
    z <- lista_srad[[i]]
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
  return(temp.dt)
  removeTmpFiles(h=0)
})

counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_rad')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
saveRDS(tabla, file=paste(counDir, '/rsrad.rds', sep=''))



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#     Load CHIRPS data          #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# R options
# Load CHIRPS data
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

###### Precipitacion
chirps_data <- list.files(pattern='*.tif$', full.names=TRUE) # lista de datos cargados con direccion
chirps_date <- list.files(pattern='*.tif$', full.names=FALSE) # quito la direccion 
chirps_date <- gsub(pattern='chirps-v2.0.', replacement='', chirps_date) # solo fecha y formato
chirps_date <- gsub(pattern='.tif', replacement='', chirps_date) #fecha
chirps_date <- strsplit(x=chirps_date, split='.', fixed=TRUE) #separo las fechas
chirps_date <- lapply(1:length(chirps_date), function(i)
{
  df <- as.data.frame(t(chirps_date[[i]]))
  colnames(df) <- c('Year', 'Month', 'Day')
  return(df)
})
chirps_date <- do.call(rbind, chirps_date)  # year, month, day
year <- as.numeric(as.character(unique(chirps_date$Year))) #year
year <- year[-length(year)]

# Create individual rasters per county

# Function to identify leap years
is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) }

# Processing CHIRPS data by county in a table ***

# Define temporal directory to save temporary rasters created during process
library(raster)
rasterOptions(tmpdir='/mnt/Workspace_cluster_9/CWR_pre-breeding/Results/Raster_t/chirps')

# Read shapefile
shp <- rgdal::readOGR(dsn = "/mnt/workspace_cluster_9/CWR_pre-breeding/Input_data/world_shape/world_shape", "all_countries")
# Calculing extet of shapefile
ext_shp <- extent(shp)

# Run process by year in parallel
library(parallel)
chirps_year <- mclapply(1:length(year), function(i){
  cat('Processing year:', year[[i]], '\n')
  # Read daily rasters in one year
  data_year <- raster::stack(chirps_data[grep(pattern=year[[i]], x=chirps_data)])
  
  # Create a rasterized version of shapefile
  template.shp <- rasterize(shp, data_year[[1]], getCover=T)
  
  # Verify if we have a leap year
  LeapYear <- is.leapyear(year[[i]])
  
  if(LeapYear==TRUE){ # Cases where we have leap years
    temp.dt <- ff(1, dim=c(ncell(data_year), 369), vmode="double")
    lapply(1:366, function(i)
    {
      z <- data_year[[i]]
      t <- getValues(z)
      t[which(t==-9999)] <- NA
      cat('Processing: Day', i, '\n')
      temp.dt[,1] <- 1:nrow(temp.dt)
      temp.dt[,2:3] <- xyFromCell(object=z, cell=1:nrow(temp.dt))
      temp.dt[,i+3] <- t[]
      return(cat("Done\n"))
    })
    temp.dt <- as.ffdf(temp.dt)
    names(temp.dt) <- c('cellID', 'lon', 'lat', paste0("d",1:366))
    temp.dt <- as.data.frame(temp.dt)
    temp.dt <- data.table(temp.dt)
    temp.dt <- temp.dt[complete.cases(temp.dt),]
  } else { # Cases where we have normal years
    temp.dt <- ff(1, dim=c(ncell(data_year), 368), vmode="double")
    lapply(1:365, function(i)
    {
      z <- data_year[[i]]
      t <- getValues(z)
      t[which(t==-9999)] <- NA
      cat('Processing: Day', i, '\n')
      temp.dt[,1] <- 1:nrow(temp.dt)
      temp.dt[,2:3] <- xyFromCell(object=z, cell=1:nrow(temp.dt))
      temp.dt[,i+3] <- t[]
      return(cat("Done\n"))
    })
    temp.dt <- as.ffdf(temp.dt)
    names(temp.dt) <- c('cellID', 'lon', 'lat', paste0("d",1:365))
    temp.dt <- as.data.frame(temp.dt)
    temp.dt <- data.table(temp.dt)
    temp.dt <- temp.dt[complete.cases(temp.dt),]
  }
  
  return(temp.dt)
}, mc.cores=20)
removeTmpFiles(h=0)

names(chirps_year) <- paste('y', year, sep='')
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/imput_tables/chirps')
if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
save(chirps_year, file=paste(counDir, '/prec_chirps.RData', sep=''))
