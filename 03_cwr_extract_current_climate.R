# CWR pre-breeding characterising testing environments: base code
# Authors: B. Mora & H. Achicanoy
# CIAT, 2017
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#                          LOAD PACKAGES
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
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
suppressMessages(if(!require(Rtsne)){install.packages('Rtsne'); library(Rtsne)} else {library(Rtsne)})
suppressMessages(if(!require(readr)){install.packages('readr'); library(readr)} else {library(readr)})
suppressMessages(if(!require(dbscan)){install.packages('dbscan'); library(dbscan)} else {library(dbscan)})
suppressMessages(if(!require(zoom)){install.packages('zoom'); library(zoom)} else {library(zoom)})

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#                       Rotar Rasters  
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

#=-=-=-=-==-=-=-=-=-=-=-=-=#
#   climate_variable_tmax    
#=-=-=-=-==-=-=-=-=-=-=-=-=#
library(parallel)
climate_variable_tmax <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmax_daily_ts_agmerra_1980_2010.nc")
climate_variable <- raster::unstack(climate_variable_tmax)
climate_variable_rotated <- mclapply(X = 1:length(climate_variable), FUN = function(i){
  rr <- raster::rotate(climate_variable[[i]])
  return(rr)
  removeTmpFiles(h = 0)
}, mc.cores = 20, mc.preschedule=F)
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax')
saveRDS(climate_variable_rotated, file=paste(counDir, '/rtmax.rds', sep=''))

#=-=-=-=-==-=-=-=-=-=-=-=-=#
#   climate_variable_tmin    
#=-=-=-=-==-=-=-=-=-=-=-=-=#
library(parallel)
climate_variable_tmin <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/tmin_daily_ts_agmerra_1980_2010.nc")
climate_variable <- raster::unstack(climate_variable_tmin)
climate_variable_rotated <- mclapply( X = 1:length(climate_variable), FUN = function(i){
  rr <- raster::rotate(climate_variable[[i]])
  return(rr)
  removeTmpFiles(h = 0)
}, mc.cores = 20, mc.preschedule=F)
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
saveRDS(climate_variable_rotated, file=paste(counDir, '/rtmin.rds', sep=''))

#=-=-=-=-==-=-=-=-=-=-=-=-=#
#   climate_variable_srad   
#=-=-=-=-==-=-=-=-=-=-=-=-=#
library(parallel)
climate_variable_srad   <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/srad_daily_ts_agmerra_1980_2010.nc")
climate_variable <- raster::unstack(climate_variable_srad  )
climate_variable_rotated <- mclapply( X = 1:length(climate_variable), FUN = function(i){
  rr <- raster::rotate(climate_variable[[i]])
  return(rr)
  removeTmpFiles(h = 0)
}, mc.cores = 20, mc.preschedule=F)
counDir <- paste('/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad')
saveRDS(climate_variable_rotated, file=paste(counDir, '/rsrad.rds', sep=''))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#                       Construcción de Tablas                                        #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#=-=-=-=#
#  TMIN #
#=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=#
#  Load Raster rotate #
#=-=-=-=-=-=-=-=-=-=-=#
tmin <-readRDS('//dapadfs/Workspace_cluster_9//CWR_pre-breeding/Results/input_tables/agmerra_tmin/rtmin.rds')
tmin <- raster::stack(tmin)
#=-=-=-=-=-=#
# To break  #
#=-=-=-=-=-=#
tmin1 <- tmin[[1:4000]]
tmin2 <- tmin[[4001:8000]]
tmin3 <- tmin[[8001:11323]]
#=-=-=-=-=-=#
#   Tables  #
#=-=-=-=-=-=#
tmin1.1 <- as.data.frame(tmin1, xy=TRUE)
tmin2.1 <- as.data.frame(tmin2, xy=TRUE)
tmin3.1 <- as.data.frame(tmin3, xy=TRUE)
#=-=-=-=-=-=#
#Save Tables#
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
saveRDS(tmin1.1, file=paste(counDir, '/tabla_tmin1.rds', sep=''))
saveRDS(tmin2.1, file=paste(counDir, '/tabla_tmin2.rds', sep=''))
saveRDS(tmin3.3, file=paste(counDir, '/tabla_tmin3.rds', sep=''))
#=-=-=-=-=-=#
#Load Tables#
#=-=-=-=-=-=#
tmin1 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tabla_tmin1.rds')
tmin2 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tabla_tmin2.rds')
tmin3 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tabla_tmin3.rds')
#=-=-=-=-=-=#
# Omit NA   #
#=-=-=-=-=-=#
system.time(tmin1.1<-tmin1[which(rowSums(is.na(tmin1))!=4000),])
system.time(tmin1.2<-tmin2[which(rowSums(is.na(tmin2))!=4000),])
#=-=-=-=-=-=#
# important #
#=-=-=-=-=-=#
rm (tmin1,tmin2)
system.time(tmin1.3<-tmin3[which(rowSums(is.na(tmin3))!=3325),])
#=-=-=-=-=-=#
# important #
#=-=-=-=-=-=#
rm(tmin3)
#=-=-=-=-=-=-=-=-=-=-=-=-#
# quit year 1980 bisiesto#
#=-=-=-=-=-=-=-=-=-=-=-=-#
tmin1.1 <- tmin1.1[,-(3:368)]

#=-=-=-=-=-=#
#Love Table # 
#=-=-=-=-=-=#
system.time(tabla <- Reduce(merge, list(tmin1.1,tmin1.2,tmin1.3)))
rm (tmin1.1,tmin1.2,tmin1.3)
#=-=-=-=-=-=#
#Save Table #
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
saveRDS(tabla, file=paste(counDir, '/table_final.rds', sep=''))
rm(tabla)

#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=#
#   TMAX    #
#=-=-=-=-=-=#

#=-=-=-=-=-=-=-=-=-=-=-=-#
#    Load Raster rotate  #
#=-=-=-=-=-=-=-=-=-=-=-=-#
tmax <-readRDS('//dapadfs/Workspace_cluster_9//CWR_pre-breeding/Results/input_tables/agmerra_tmax/rtmax.rds')
tmax <- raster::stack(tmax)
#=-=-=-=-=-=#
# To break
#=-=-=-=-=-=#
tmax1 <- tmax[[1:4000]]
tmax2 <- tmax[[4001:8000]]
tmax3 <- tmax[[8001:11323]]

#=-=-=-=-=-=#
#  Tables
#=-=-=-=-=-=#
tmax1.1 <- as.data.frame(tmax1, xy=TRUE)
tmax2.1 <- as.data.frame(tmax2, xy=TRUE)
tmax3.1 <- as.data.frame(tmax3, xy=TRUE)

#=-=-=-=-=-=#
#Save Tables
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax')
saveRDS(tmax1.1, file=paste(counDir, '/tabla_tmax1.rds', sep=''))
saveRDS(tmax2.1, file=paste(counDir, '/tabla_tmax2.rds', sep=''))
saveRDS(tmax3.1, file=paste(counDir, '/tabla_tmax3.rds', sep=''))
rm(tmax)
#=-=-=-=-=-=#
#Load Tables
#=-=-=-=-=-=#
tmax1 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tabla_tmax1.rds')
tmax2 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tabla_tmax2.rds')
tmax3 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tabla_tmax3.rds')

#=-=-=-=-=-=#
#  Omit NA
#=-=-=-=-=-=#
system.time(tmax1.1<-tmax1[which(rowSums(is.na(tmax1))!=4000),])
system.time(tmax1.2<-tmax2[which(rowSums(is.na(tmax2))!=4000),])
#=-=-=-=-=-=#
# important
#=-=-=-=-=-=#
rm (tmax1,tmax2)
system.time(tmax1.3<-tmax3[which(rowSums(is.na(tmax3))!=3325),])
#=-=-=-=-=-=#
# important
#=-=-=-=-=-=#
rm(tmax3)
#=-=-=-=-=-==-=-=-=-=--=-#
#quit year 1980 bisiesto
#=-=-=-=-=-==-=-=-=-=--=-#
tmax1.1 <- tmax1.1[,-(3:368)]
#=-=-=-=-=-=#
# Love Table 
#=-=-=-=-=-=#
system.time(tabla <- Reduce(merge, list(tmax1.1,tmax1.2,tmax1.3)))
rm (tmax1.1,tmax1.2,tmax1.3)
#=-=-=-=-=-=#
#Save Table
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin')
saveRDS(tabla, file=paste(counDir, '/table_final.rds', sep=''))
rm(tabla)
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=#
#   SRAD
#=-=-=-=-=-=#

#=-=-=-=-=-=-=-=-=-=-=-=-#
#   Load Raster rotate   #
#=-=-=-=-=-=-=-=-=-=-=-=-#
srad <-readRDS('//dapadfs/Workspace_cluster_9//CWR_pre-breeding/Results/input_tables/agmerra_srad/rsrad.rds')
srad <- raster::stack(srad)
#=-=-=-=-=-=#
#  To break
#=-=-=-=-=-=#
srad1 <- srad[[1:4000]]
srad2 <- srad[[4001:8000]]
srad3 <- srad[[8001:11323]]

#=-=-=-=-=-=#
#   Tables 
#=-=-=-=-=-=#
srad1.1 <- as.data.frame(srad1, xy=TRUE)
srad2.1 <- as.data.frame(srad2, xy=TRUE)
srad3.1 <- as.data.frame(srad3, xy=TRUE)

#=-=-=-=-=-=#
#Save Tables#
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad')
saveRDS(srad1.1, file=paste(counDir, '/tabla_srad1.rds', sep=''))
saveRDS(srad2.1, file=paste(counDir, '/tabla_srad2.rds', sep=''))
saveRDS(srad3.1, file=paste(counDir, '/tabla_srad3.rds', sep=''))
rm(srad)

#=-=-=-=-=-=#
#Load Tables
#=-=-=-=-=-=#
srad1 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/tabla_srad1.rds')
srad2 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/tabla_srad2.rds')
srad3 <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/tabla_srad3.rds')

#=-=-=-=-=-=#
#  Omit NA
#=-=-=-=-=-=#
system.time(srad1.1<-srad1[which(rowSums(is.na(srad1.1))!=4000),])
system.time(srad2.1<-srad2[which(rowSums(is.na(srad2.1))!=4000),])
#=-=-=-=-=-=#
# important
#=-=-=-=-=-=#
rm (srad1.1,srad2.1)
system.time(srad3.1<-srad3[which(rowSums(is.na(srad3))!=3323),])
#=-=-=-=-=-=#
# important
#=-=-=-=-=-=#
rm(srad3)
#=-=-=-=-=-==-=-=-=-=--=-#
# quit year 1980 bisiesto
#=-=-=-=-=-==-=-=-=-=--=-#
srad1.1 <- srad1.1[,-(3:368)]
#=-=-=-=-=-=#
# Love Table 
#=-=-=-=-=-=#
system.time(tabla <- Reduce(merge, list(srad1.1,srad2.1,srad3.1)))

rm (srad1.1,srad2.1,srad3.1)
#=-=-=-=-=-=#
# Save Table
#=-=-=-=-=-=#
counDir <- paste('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad')
saveRDS(tabla, file=paste(counDir, '/table_final.rds', sep=''))
rm(tabla)

#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#     Load CHIRPS data          #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)
chirps_dir <- Sys.info(); chirps_dir <- chirps_dir[names(chirps_dir)=="sysname"]
if(chirps_dir == "Linux"){
  wk_dir <- "/mnt/data_cluster_4/observed/gridded_products/chirps/daily/32bits"; setwd(wk_dir); rm(wk_dir)
} else {
  if(chirps_dir == "Windows"){
    wk_dir <- "//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily/32bits"; setwd(wk_dir); rm(wk_dir)
  }
}; 
rm(chirps_dir)

library(parallel)
##Loop with lapply everybody  years##
year <- 1981:2010
chirps_data_prec <- list.files(pattern=".tif$", full.names=TRUE)
system.time(data_day <- mclapply(1:length(year), function(i){ 
  cat(paste("Procesando year:", year[i], "\n", sep = "")) 
  rasters  <- raster::stack(chirps_data_prec[grep(pattern=year[[i]], x=chirps_data_prec)])
  return(rasters)
  removeTmpFiles(h=0)
}, mc.cores = 10))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#            PRUEBA             #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

for(i in 1:365){ clase<-class(pre_year24[[i]]);cat("object",clase," ",i,"\n");Sys.sleep(0.2)}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#proceso debe hacerse para cada año por separado# 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
year1<- data_day[[1]]
year1<- raster::unstack(year1)
pre_year1<- mclapply(X = 1:length(year1), FUN = function(i){
  rr <- raster::resample(year1[[i]],modo)
  return(rr)
  removeTmpFiles(h = 0)
}, mc.cores = 10, mc.preschedule=F)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#            PRUEBA             #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

for(i in 1:365){ clase<-class(pre_year1[[i]]);cat("object",clase," ",i,"\n");Sys.sleep(0.2)}
year1<- raster::stack(pre_year1)
tabla<- as.data.frame(year1, xy=TRUE)
system.time(tabla1<-tabla[which(rowSums(is.na(tabla))!=365),])
saveRDS(tabla1,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/tabla1.rds")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#Despues de hacer esto para cada año por separado #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=-=-=-=-=-=-#
#      quito los -9999
#=-=-=-=-=-=-=-=-=-=-=-#

g<- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)
chirps_dir <- Sys.info(); chirps_dir <- chirps_dir[names(chirps_dir)=="sysname"]
if(chirps_dir == "Linux"){
  wk_dir <- "/mnt/workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps"; setwd(wk_dir); rm(wk_dir)
} else {
  if(chirps_dir == "Windows"){
    wk_dir <- "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps"; setwd(wk_dir); rm(wk_dir)
  }
}; 
rm(chirps_dir)
chirps_data_prec <- list.files(pattern=".rds$", full.names=TRUE)

library(gtools)
chirps_data <- mixedsort(chirps_data_prec)






tabla <- 1:30
tab <- lapply(1:length(tabla), function(i) {  
  cat(paste("Procesando year:", tabla[i], "\n", sep = "")) 
  t <-  readRDS(chirps_data_prec[i])
  
  if (ncol(t)==367) {
    (tabla<-t[which(rowSums(t) > (0)),])
  } 
  
  if (ncol(t) == 368)   
  {
    (tabla<-t[which(rowSums(t) > (0)),])
  }
  
  return(t)
  removeTmpFiles(h = 0)
})

t1 <- tab[[1]]
t2 <- tab[[2]]
t3 <- tab[[3]]
t4 <- tab[[4]]
t5 <- tab[[5]]
t6 <- tab[[6]]
t7 <- tab[[7]]
t8 <- tab[[8]]
t9 <- tab[[9]]
t10 <- tab[[10]]
t11 <- tab[[11]]
t12 <- tab[[12]]
t13 <- tab[[13]]
t14 <- tab[[14]]
t15 <- tab[[15]]

system.time(tabla1.1 <- Reduce(merge, list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15)))
saveRDS(tabla1.1, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/tabla1.1.rds")

#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#

t16 <- tab[[16]]
t17 <- tab[[17]]
t18 <- tab[[18]]
t19 <- tab[[19]]
t20 <- tab[[20]]
t21 <- tab[[21]]
t22 <- tab[[22]]
t23 <- tab[[23]]
t24 <- tab[[24]]
t25 <- tab[[25]]
t26 <- tab[[26]]
t27 <- tab[[27]]
t28 <- tab[[28]]
t29 <- tab[[29]]
t30 <- tab[[30]]

system.time(tabla1.2<- Reduce(merge, list(t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30)))
saveRDS(tabla1.2, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/tabla1.2.rds")

tabla_final <- Reduce(merge, list(tabla1.1, tabla1.2))
saveRDS(tabla_final, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/tabla_final.rds")

#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-==--=-=-==-==-=-=-=-===-=-=-=-=-=-=-=-=-=-=-#






