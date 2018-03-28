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
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(FactoClass)){install.packages('FactoClass'); library(FactoClass)} else {library(FactoClass)})
suppressMessages(if(!require(ade4)){install.packages('ade4'); library(ade4)} else {library(ade4)})
suppressMessages(if(!require(xtable)){install.packages('xtable'); library(xtable)} else {library(xtable)})
suppressMessages(if(!require(ggdendro)){install.packages('ggdendro'); library(ggdendro)} else {library(ggdendro)})
suppressMessages(if(!require(compiler)){install.packages('compiler'); library(compiler)} else {library(compiler)})
suppressMessages(if(!require(ggthemes)){install.packages('ggthemes'); library(ggthemes)} else {library(ggthemes)})
#suppressMessages(if(!require(dtwclust)){install.packages('dtwclust'); library(dtwclust)} else {library(dtwclust)})
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
if(OSys == "Linux"){
  root <- "/mnt/workspace_cluster_9"
  base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
} else {
  if(OSys == "Windows"){
    root <- "//dapadfs/Workspace_cluster_9"
    base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
  }
}; rm(OSys)




# Brayan Mora
future <- list.files(path = "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k", pattern = "*.nc$", full.names = F)
gcm1_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp45_r1i1p1_CCSM4", x = future)]
gcm2_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp45_r1i1p1_MIROC-ESM-CHEM", x = future)]
gcm3_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_ACCESS1-0", x = future)]
gcm4_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_GFDL-ESM2M", x = future)]
gcm5_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_inmcm4.", x = future)]

g <- gcm1_prec[-(1:14)]
g1 <- g[-(51:length(g))]


library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

cores<- detectCores()
cl<- makeCluster(cores-20)
registerDoParallel(cl) 


# mylist <- rep(list(NA),50)

system.time(l <- foreach(i=1:50) %dopar% {
  mylist <- list()
  root <- "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k/"
  require(dplyr)
  require(velox)
  require(raster)
  require(raster)
  
  cat(paste0("procesando tabla:",i,"\n" ))
  prec <- raster::stack(paste0(root, g[i]))  
  prec_vx <- velox::velox(prec)
  pnt <- prec_vx$getCoordinates()
  colnames(pnt) <- c("x", "y")
  pnt <- data.frame(pnt)
  pnt$x[pnt$x > 180] <- pnt$x[pnt$x > 180] - 360
  rtmp <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/AgMerra_template.RDS")
  rtmp[which(!is.na(rtmp[]))] <- 1
  # rtmp[which(!is.na(rtmp[]))] <- 1
  # pnt_rtmp <- coordinates(rtmp[which(!is.na(rtmp[]))])
  pnt_rtmp <- rasterToPoints(rtmp)
  pnt_rtmp <- pnt_rtmp[,1:2]
  pnt_rtmp <- data.frame(pnt_rtmp)
  pnt_rtmp$cellID <- cellFromXY(object = rtmp, xy = pnt_rtmp)
  pnt$cellID <- cellFromXY(object = rtmp, xy = pnt)
  pnt <- pnt[pnt$cellID %in% pnt_rtmp$cellID,]
  rownames(pnt) <- 1:nrow(pnt)
  pnt$x[pnt$x <= 0] <- pnt$x[pnt$x <= 0] + 360
  vls <- prec_vx$extract_points(sp = SpatialPoints(pnt[,1:2]))
  vls <- as_data_frame(vls)
  vls <- cbind(pnt, vls)
  vls$x[vls$x > 180] <- vls$x[vls$x > 180] - 360
  vls1<- data.frame(cellID= vls$cellID, lon= vls$x, lat= vls$y, vls[,-3])
  vls1$x <- NULL 
  vls1$y <- NULL
  names(vls1)[4:ncol(vls1)] <- names(prec)
  names(vls1)[4:ncol(vls1)] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(vls1)[4:ncol(vls1)]), fixed = T))
  mylist[[i]] <- vls1
  
  
} )
stopCluster(cl)

l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:3)]
  return(x)
})
tmin  <- do.call(cbind,l1)
tmin <- cbind(l[[1]], tmin)


##Filtro

america <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_america.rds")
africa  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_africa.rds")
oceania <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_oceania.rds")
asia    <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_asia.rds")
europa  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_europa.rds")



a <- c(america,africa,oceania,asia,europa)
filtro <- lapply(1:length(a), function(i){
  tabla <- dplyr::filter(tmin, cellID %in% a[i]$cellID)
  tabla[]
  return(tabla)
})

cell <- ch$cellID

tminf<- dplyr::filter(tmin, cellID %in% cellreferencia)



### conversion  

system.time(saveRDS(tmin, paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_future_climate/rcp45/gcm5/chirps/pr_day_BCSD_rcp85_r1i1p1_inmcm4.rds")))
