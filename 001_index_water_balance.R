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


rbalance<-  function(continent= "Africa"){
  
  output <- paste0(root, "/CWR_pre-breeding/Results/Potato/Hydric_balance/potato_water_indices_", tolower(continent), ".rds")
  if(!file.exists(output)){
    prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_', tolower(continent), '.rds'))
    crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Potato/database/area_base.rds"))
    require(dplyr)
    prec_filter <- dplyr::filter(crop_area, cellID %in% prec$cellID)
    prec_filter <- prec_filter[!duplicated(prec_filter$cellID),]
    cpc <- unique(prec_filter$cellID)
    
    
    
    system.time(wb <- lapply(1 :length(cpc), function(i){# 
      require( tidyr)
      cat(paste0("looking for coordinate: ", i, "\n"))
      wb <- list.files(path = paste0(root, "/CWR_pre-breeding/Input_data/_soils/Water_balance/",tolower(continent)) , paste0(cpc[i],".rds"), full.names = T)
      wb1 <- lapply(wb, readRDS)
      wb1<- do.call(rbind,wb1)
      wb1 <- data.frame(Date = row.names(wb1), cellID = wb1$cellID, lon= wb1$lon, lat= wb1$lat, ERATIO= wb1$ERATIO)
      row.names(wb1)<- NULL
      
      if(dim(wb1)==0){
        cat(paste0("archivo malo for coordinate: ", i, "\n"))
        return(wb1)
      }else{
        wb1 <- wb1 %>% spread("Date","ERATIO")
        return(wb1)}
    }))
    
    wb <- do.call(rbind, wb)
    cat(">>> Loading crop cycle data ...\n")
    
    # Planting dates
    planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
    planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
    # Harversting dates
    harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
    harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
    
    
    # Extract important dates
    wb$Planting <- raster::extract(x = planting_rf_ggcmi, y = wb[,c("lon", "lat")]);
    wb$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = wb[,c("lon", "lat")]); 
    wb$Duration <- ifelse(test = wb$Planting < wb$Harvest, yes = "One year", no = "Two years")
    wb$Harvest[which(wb$Harvest == -99)] <- NA
    wb<- na.omit(wb)
    wb$Harvest <- round(wb$Harvest)
    wb$Planting <- round(wb$Planting)
    
    test<- wb[which(wb$Duration =="Two years"), ]
    
    
    
    library(doSNOW)
    library(foreach)
    library(parallel)
    library(doParallel)
    
    cores<- detectCores()
    cl<- makeCluster(cores-18)
    registerDoParallel(cl) 
    
    system.time(l <- foreach(i=1:nrow(wb)) %dopar% {
      mylist <- list()
      duration <- wb$Duration[i]
      start <- wb$Planting[i]
      end <- wb$Harvest[i]
      time.serie <- wb[i, 1:(ncol(wb)-3)]
      if(duration == "One year"){
        require(dplyr)
        require(tidyr)
        require(compiler)
        
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start +70& Yday <= end)
        
        # Drought
        dr_stress<- function(wb, p_thresh = 0.5){
          days <- sum(wb < 0.5)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        drought <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(drought = dr_stressCMP(Value, p_thresh= 0.5))
        drought <- drought %>% as.data.frame
        names(drought)[2] <- "Value"; drought$Variable <- "drought"
        drought <- data.frame(cellID= unique(X$cellID), drought)
        
        
      }else{
        require(dplyr)
        require(tidyr)
        require(compiler)
        
        
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% filter(Yday %in% c(258:365, 1:end))
        X <- X[-(1:(end)),]
        X <- X[-((nrow(X)-(365-start)): nrow(X)),]
        X$condition <- NA
        try(for(j in 1:nrow(X)){
          X$condition[j] <- X$Yday[j+1] - X$Yday[j]
        })
        c <- c(unique(X$condition))
        X$condition[which(is.na(X$condition))] <- c[3]
        chngs <- which(X$condition == c[3])
        for(j in 1:length(chngs)){
          if(j == 1){
            X$condition2[1:chngs[j]] <- j
          } else {
            X$condition2[(chngs[j-1]+1):(chngs[j])] <- j
          }
        }
        X$condition2[is.na(X$condition2)] <- length(chngs)
        X$condition <- NULL
        names(X)[ncol(X)] <- "condition"
        rownames(X)<- 1:nrow(X)
        
        corte <- c(1:29)
        lista <- lapply(1:length(corte), function(i){
          tablas <- filter (X, condition== corte[i]) 
          return(tablas)  
        })
        
        tab <- lapply(1:length(lista), function(i){
          tabla <-lista[[i]][(70:nrow(lista[[i]])),]
          return(tabla)
        })
        X <- do.call(rbind,tab)
        X <- na.omit(X)
        
        
        # Drought
        dr_stress<- function(wb, p_thresh = 0.5){
          days <- sum(wb < 0.5)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        drought <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(drought = dr_stressCMP(Value, p_thresh= 0.5))
        drought <- drought %>% as.data.frame
        names(drought)[2] <- "Value"; drought$Variable <- "drought"
        drought <- data.frame(cellID= unique(X$cellID), drought)
        
        
      }
      
      mylist[[i]] <- drought
      
    } )
    stopCluster(cl)
    
    tabla <- do.call(rbind, l)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    return(cat("Process done\n"))
    
  }
  else{
    cat(">>> Indices have been already calculated ...\n")
  }
}

system.time( rbalance(continent = "America"))
