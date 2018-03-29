#CWR pre-breeding characterising testing environments: Index sunflower 
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

# Load packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(FactoClass)){install.packages('FactoClass'); library(FactoClass)} else {library(FactoClass)})
suppressMessages(if(!require(ade4)){install.packages('ade4'); library(ade4)} else {library(ade4)})
suppressMessages(if(!require(xtable)){install.packages('xtable'); library(xtable)} else {library(xtable)})
suppressMessages(if(!require(ggdendro)){install.packages('ggdendro'); library(ggdendro)} else {library(ggdendro)})
suppressMessages(if(!require(compiler)){install.packages('compiler'); library(compiler)} else {library(compiler)})
suppressMessages(if(!require(ggthemes)){install.packages('ggthemes'); library(ggthemes)} else {library(ggthemes)})
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
suppressMessages(if(!require(googlesheets)){install.packages("googlesheets");library(googlesheets)}else{library(googlesheets)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})
suppressMessages(if(!require(caTools)){install.packages("caTools");library(caTools)}else{library(caTools)})
library(parallel)

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
}; rm(OSys)



waterbalance<-  function(continent= "Oceania", ncores= 10 ){
  
  output <- paste0(root, "/CWR_pre-breeding/Results/Potato/Hydric_balance/potato_water_indices_", tolower(continent), ".rds")
  if(!file.exists(output)){
    prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_', tolower(continent), '.rds'))
    crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Potato/database/area_base.rds"))
    require(dplyr)
    prec_filter <- dplyr::filter(crop_area, cellID %in% prec$cellID)
    prec_filter <- prec_filter[!duplicated(prec_filter$cellID),]
    cpc <- unique(prec_filter$cellID)
    
    
    wb <- mclapply(1 : length(cpc), function(i){
      cat(paste0("looking for coordinate: ", i, "\n"))
      wb <- list.files(path = paste0(root, "/CWR_pre-breeding/Input_data/_soils/Water_balance/",tolower(continent)) , paste0(cpc[i],".rds"), full.names = T)
      wb1 <- lapply(wb, readRDS)
      wb1<- do.call(rbind,wb1)
      return(wb1)
    },mc.cores = ncores, mc.preschedule = F)
    
    
    cat(">>> Loading crop cycle data ...\n")
    
    # Planting dates
    planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
    planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
    # Harversting dates
    harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
    harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
    
    
    cat(">>> Calculating indices ...\n")
    
    index <-  mclapply(1:length(wb), function(j){
      iwb <- wb[[j]]
      iwb$SRAD <- NULL
      iwb$TMIN <- NULL
      iwb$TMAX <- NULL
      iwb$RAIN <- NULL
      iwb$CUM_RAIN <- NULL
      iwb$DEMAND <- NULL
      iwb$RUNOFF <- NULL
      iwb$AVAIL <- NULL
      iwb$ETMAX <- NULL
      iwb$Planting <- raster::extract(x = planting_rf_ggcmi, y = iwb[,c("lon", "lat")])
      iwb$Harvest <- raster::extract(x = harvest_rf_ggcmi, y = iwb[,c("lon", "lat")])
      iwb$Duration <- ifelse(test =iwb$Planting < iwb$Harvest , yes = "One year", no = "Two years")
      
      duration <- unique(iwb$Duration)
      start <-  unique(iwb$Planting)
      end <- unique(iwb$Harvest)
      
      
      if(duration == "One year"){
        
        iwb$Date<- row.names(iwb)
        row.names(iwb)<- NULL
        
        iwb$condition <- NA
        iwb$condition[which(iwb$ERATIO < 0.5 )] <- 1
        iwb$condition[which(iwb$ERATIO >= 0.5 )] <- 0
        
        
        iwb$Year <- lubridate::year(as.Date(iwb$Date))
        iwb$Yday <- lubridate::yday(as.Date(iwb$Date))
        iwb <- iwb%>% group_by(Year) %>% dplyr::filter(Yday >= start+70 & Yday <= end)
        
        
        # Drought
        dr_stress<- function(condition, p_thresh = 1){
          days <- sum(condition >= 1)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        drought <- iwb %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(drought = dr_stressCMP(condition, p_thresh= 1))
        
        drought <- drought %>% as.data.frame
        names(drought)[2] <- "Value"; drought$Variable <- "drought"
        
        drought <- data.frame(cellID= unique(iwb$cellID), drought)
        return(drought)
        
      }else{
        
        iwb$Date<- row.names(iwb)
        row.names(iwb)<- NULL
        iwb$con<- NA
        iwb$condition[which(iwb$ERATIO < 0.5 )] <- 1
        iwb$condition[which(iwb$ERATIO >= 0.5 )] <- 0
        
        
        
        iwb$Year <- lubridate::year(as.Date(iwb$Date))
        iwb$Yday <- lubridate::yday(as.Date(iwb$Date))
        
        iwb <- iwb %>% filter(Yday %in% c(start:365, 1:end))
        iwb <- iwb[-(1:(end)),]
        iwb <- iwb[-((nrow(iwb)-(365-start)): nrow(iwb)),]
        iwb$condition <- NA
        for(j in 1:nrow(iwb)){
          iwb$condition[j] <- iwb$Yday[j+1] - iwb$Yday[j]
        }
        c <- c(unique(iwb$condition))
        iwb$condition[which(is.na(iwb$condition))] <- c[3]
        chngs <- which(iwb$condition == c[3])
        for(j in 1:length(chngs)){
          if(j == 1){
            iwb$condition2[1:chngs[j]] <- j
          } else {
            iwb$condition2[(chngs[j-1]+1):(chngs[j])] <- j
          }
        }
        iwb$condition2[is.na(iwb$condition2)] <- length(chngs)
        iwb$condition <- NULL
        names(iwb)[ncol(iwb)] <- "condition"
        rownames(iwb)<- 1:nrow(iwb)
        
        corte <- c(1:29)
        lista <- lapply(1:length(corte), function(i){
          tablas <- filter (X, condition== corte[i]) 
          return(tablas)  
        })
        
        tab <- lapply(1:length(lista), function(i){
          tabla <-lista[[i]][(70:end),]
          return(tabla)
        })
        iwb <- do.call(rbind,tab)
        iwb <- na.omit(iwb)
        
        # Drought
        dr_stress<- function(condition, p_thresh = 1){
          days <- sum(condition >= 1)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        drought <- iwb %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(drought = dr_stressCMP(condition, p_thresh= 1))
        
        drought <- drought %>% as.data.frame
        names(drought)[2] <- "Value"; drought$Variable <- "drought"
        
        drought <- data.frame(cellID= unique(iwb$cellID), drought)
        return(drought)
      }
      
    },mc.cores = ncores, mc.preschedule = F)
    
    tabla <- do.call(rbind, index)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    return(cat("Process done\n"))
    
  }
  else{
    cat(">>> Indices have been already calculated ...\n")
  }
}
