# Validating waterlogging stress indices - CWR pre-breeding
# Author: H. Achicanoy
# CIAT, 2018

# Load libraries
# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)

# Load packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(tcltk)){install.packages("tcltk");library(tcltk)}else{library(tcltk)})

# Path settings
OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
if(OSys == "Linux"){
  root <- "/mnt/workspace_cluster_9"
  base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
} else {
  if(OSys == "Windows"){
    root <- "//dapadfs/Workspace_cluster_9"
    base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
  }
}

# Load drought stress indices by crop-stress and continent
# (It should vary by crop, stress and continent)

# seleccionar continente
continent <- "Asia"
gcm <- "gcm5"
#######################################################
#######################################################
#######################################################

prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_', tolower(continent), '.rds'))
crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Bean/database/area_base.rds"))
require(dplyr)
prec_filter <- dplyr::filter(crop_area, cellID %in% prec$cellID)
prec_filter <- prec_filter[!duplicated(prec_filter$cellID),]
px_list <- c (unique(prec_filter$cellID))


# Load crop cycle calendar by crop
# Planting date
planting_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harvest day
harvest_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

# Load historical water balance for pixels where indices have been calculated


wb_files <- paste0(root, "/CWR_pre-breeding/Input_data/_soils/Water_balance_future/rcp85/",gcm,"/",continent,"/cellID_",px_list,".rds")
wb_files <- wb_files[which(file.exists(wb_files))]
px_list  <- px_list[which(file.exists(wb_files))]


# Calculating number of days within crop cycle with ERATIO < 0.5 (evidence of drought)

library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

cores<- detectCores()
cl<- makeSOCKcluster(cores-25)
registerDoParallel(cl) 


ntasks <- length(px_list)
pb <- tkProgressBar(max=ntasks)
progress <-  function(n)
  setTkProgressBar(pb, n)

opts <- list(progress=progress)



require(parallel)

system.time(hidric <- foreach(i=1:length(px_list), .options.snow=opts,.packages='dplyr') %dopar% { 
  tryCatch({ mylist <- list()
  wb <- readRDS(wb_files[i])
  wb$watbal[[1]]$Year <- lubridate::year(as.Date(rownames(wb$watbal[[1]])))
  wb$watbal[[1]]$Yday <- lubridate::yday(as.Date(rownames(wb$watbal[[1]])))
  
  # Extracting crop cycle
  start <- raster::extract(x = planting_rf_ggcmi, y = wb$watbal[[1]][,c("lon", "lat")] %>% unique)
  end   <- raster::extract(x = harvest_rf_ggcmi, y = wb$watbal[[1]][,c("lon", "lat")] %>% unique)
  duration <- ifelse(test = start< end, yes = "One year", no = "Two years")
  
  
  if(duration== "One year"){
    require(dplyr)
    runoff <- wb$watbal[[1]] %>%
      group_by(Year) %>%
      dplyr::filter(Yday >= start & Yday <= end-2) %>%
      dplyr::arrange(Yday) %>%
      dplyr::summarise(Runoff = sum(RUNOFF > 0, na.rm = T)) %>%
      as.data.frame
    runoff$cellID <- unique(as.numeric(wb$cellID))
    
    
  }else{
    require(dplyr)
    wb  <- wb$watbal[[1]]
    wb <- wb %>% filter(Yday %in% c(start:365, 1:end))
    wb <- wb[-(1:(which(wb$Yday == end-1)[1])),]
    wb <- wb[-(which(wb$Yday == start)[length(which(wb$Yday == start))]:nrow(wb)),]
    
    wb$condition <- NA
    try(for(j in 1:nrow(wb)){
      wb$condition[j] <- wb$Yday[j+1] - wb$Yday[j]
    })
    c <- c(unique(wb$condition))
    c <- sort(c) 
    wb$condition[which(is.na(wb$condition))] <- c[4]
    chngs <- which(wb$condition == c[4])
    
    row.names(wb) <- 1:nrow(wb)
    # chngs <- chngs[-1]
    for(j in 1:length(chngs)){
      if(j == 1){
        wb$condition2[1:chngs[j]] <- j
      } else {
        wb$condition2[(chngs[j-1]+1):(chngs[j])] <- j
      }
    }
    
    wb$condition2[is.na(wb$condition2)] <- length(chngs)
    wb$condition <- NULL
    names(wb)[ncol(wb)] <- "condition"
    rownames(wb)<- 1:nrow(wb)
    
    runoff <- wb %>% group_by(condition) %>%
      dplyr::summarise(Runoff = sum(RUNOFF > 0, na.rm = T))%>%
      as.data.frame
    lon <- nrow(runoff)
    if (lon== 30){
      runoff <- data.frame(cellID = unique(as.numeric(wb$cellID)), Year = unique(wb$Year),Runoff = runoff$Runoff )
    }else{
      runoff <- data.frame(cellID = unique(as.numeric(wb$cellID)), Year = unique(wb$Year)[-1],Runoff = runoff$Runoff )
    }
  }
  mylist[[i]] <- runoff
  }, error=function(e){})
  
})
stopCluster(cl)

results <- do.call(rbind, hidric)
head(results); tail(results)
saveRDS(results,paste0 (root,"/CWR_pre-breeding/Results/Bean/Waterlogging/_future/rcp85/",gcm,"/waterlogging_rbind_",continent,".rds" ))
