# CWR pre-breeding characterising testing environments: Index calculation for current data
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017

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

l1 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl1_10km_ll.tif"))
l2 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl2_10km_ll.tif"))
l3 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl3_10km_ll.tif"))
l4 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl4_10km_ll.tif"))
l5 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl5_10km_ll.tif"))
l6 <- raster::stack(paste0(root,"/CWR_pre-breeding/Input_data/_soils/AWCh2/AWCh2_M_sl6_10km_ll.tif"))
l <- raster::stack(l1,l2,l3,l4,l5,l6)
rm(l1,l2,l3,l4,l5,l6)
r <- resample(l,base, method="ngb" )
rm(l)

##Data Presence
###############test

cor <- readRDS(paste0(root,("/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds")))


######################################################
tabla <- data.frame( raster::extract(r,cor[,c("lon","lat")]))
tabla <- data.frame(cor, tabla)
tabla$raiz <- 20

tabla <- tabla %>% select(lon, lat, raiz, 4:9)
depths <- c(25,100,225,450,800,1500)
colnames(tabla)[4:ncol(tabla)] <- paste0("d.", depths)

soilcap_calc_mod <- function(x, minval, maxval) {
  if(!is.na(x[3])){
    rdepth <- max(c(x[3],minval)) #cross check
    rdepth <- min(c(rdepth,maxval)) #cross-check
    wc_df <- data.frame(depth=c(2.5,10,22.5,45,80,150),wc=(x[4:9])*.01)
    if (!rdepth %in% wc_df$depth) {
      wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
      wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
      y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
      x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
      ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
      wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
    }
    wc_df <- wc_df[which(wc_df$depth <= rdepth),]
    wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
    wc_df$soilcap <- wc_df$soilthick * wc_df$wc
    soilcp <- sum(wc_df$soilcap) * 10 #in mm
    return(soilcp)
  } else {
    soilcp <- NA
    return(soilcp)
  }
}
tabla$soilcp <- apply(tabla, 1, FUN=soilcap_calc_mod, minval=45, maxval=100)
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/_soils/solil_capacity.rds"))
capacidad  <- data.frame(cellID= cor$cellID, lon=tabla$lon, lat = tabla$lat, cap = tabla$soilcp)
rm(tabla)
######################################################################################################
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_africa.rds'))
prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_europa.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_asia.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_oceania.rds'))


#############################################
##### DATOS DE CICLO DE CULTIVO GGCMI   #####
#############################################

# Planting dates
planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(planting_rf_ggcmi)
prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(harvest_rf_ggcmi)
prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")

# Just for verifying
prec[1:5,(ncol(prec)-5):ncol(prec)]
table(prec$Duration)
cell <- prec$cellID
cap_filter<- dplyr::filter(capacidad, capacidad$cellID %in% cell)
cap_filter <- cap_filter[!duplicated(cap_filter$cellID),]
prec1 <- dplyr::filter(prec, prec$cellID %in% cap_filter$cellID)
rm(prec,capacidad, cap_filter1)

prec1$soil_cap <- cap_filter$cap
prec1[1:5,(ncol(prec1)-5):ncol(prec1)]
prec <- prec1
rm(prec1)
## INDEX DE ANEGAMIENTO 
prec[1:5,(ncol(prec)-5):ncol(prec)]

library(parallel)
system.time( inundacion_index <- mclapply(1:nrow(prec), function(i){
  cat(paste0("Processed pixel:", i, "\n"))
  # Parameters
  duration <- prec$Duration[i]
  start <- prec$Planting[i]
  end <- prec$Harvest[i]
  
  # Just one pixel
  time.serie <- prec[i, 1:(ncol(prec)-5)]
  if(duration == "One year"){
    
    suppressMessages(library(tidyverse))
    suppressMessages(library(compiler))
    
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
    vec <- cap_filter[,"cap"]
    X$cap <- vec[i]
    X$diferencia <- (X$Value- X$cap)
    X$si <- NA
    X$si[which(X$diferencia < 0)]<- 0
    X$si[which(X$diferencia >= 0)]<- 1
    
    # Anegamiento 
    inundacion <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(inundacion = sum(si))
    inundacion <- inundacion %>% as.data.frame
    names(inundacion)[2] <- "Value"; inundacion$Variable <- "INUNDACION"
    results <- data.frame(cellID = unique(X$cellID), inundacion )
  } else {
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
    vec <- cap_filter[,"cap"]
    X$cap <- vec[i]
    X$diferencia <- (X$Value- X$cap)
    X$si <- NA
    X$si[which(X$diferencia < 0)]<- 0
    X$si[which(X$diferencia >= 0)]<- 1
    
    # Anegamiento 
    inundacion <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(inundacion = sum(si))
    inundacion <- inundacion %>% as.data.frame
    names(inundacion)[2] <- "Value"; inundacion$Variable <- "INUNDACION"
    
    results <- data.frame(cellID = unique(X$cellID), inundacion )
    
  }
  return(results)
}, mc.cores = 10, mc.preschedule = F))
tabla <- do.call(rbind,inundacion_index)
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/inundacion_america.rds"))
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/inundacion_africa.rds"))
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/inundacion_asia.rds"))
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/inundacion_europa.rds"))
saveRDS(tabla, paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/inundacion_oceaniards"))




