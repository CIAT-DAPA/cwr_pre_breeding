# Validating waterlogging stress indices - CWR pre-breeding
# Author: H. Achicanoy
# CIAT, 2018

# Load libraries
# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)


start <- 1 
end <- 90 



calculate_vpd <- function(continent ="Oceania", start= 360, end = 30, ncores= 15){
  
  
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
  
  output <- paste0(root, "/CWR_pre-breeding/Proyect_bean/Calculate_vpd/_current/median_vpd_", tolower(continent), ".rds")
  if(!file.exists(output)){
    
    #Para_filtrar_chirps_ggcmi
    prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_', tolower(continent), '.rds'))
    crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Bean/database/area_base.rds"))
    require(dplyr)
    prec_filter <- dplyr::filter(crop_area, cellID %in% prec$cellID)
    prec_f <- prec_filter[!duplicated(prec_filter$cellID),]
    rm(prec, crop_area,prec_filter)
    
    # Cargar_entradas
    tmin <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_",continent,".rds"))
    tmax <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmax/tmax_filtered_",continent,".rds"))
    srad <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/_current_climate/agmerra_srad/srad_filtered_",continent,".rds"))
    #Filtrar_entradas
    tmin <- dplyr::filter(tmin, cellID %in% prec_f$cellID)
    tmax <- dplyr::filter(tmax, cellID %in% prec_f$cellID)
    srad <- dplyr::filter(srad, cellID %in% prec_f$cellID)
    
    
    require(parallel)
    system.time(indexes_been <- mclapply(1:nrow(tmax), function(i){   ### nrow(tmax)
      
      time.serie  <- tmax[i, 1:(ncol(tmax)-1)]
      time.serie1 <- tmin[i, 1:(ncol(tmin)-1)]
      time.serie2 <- srad[i, 1:(ncol(srad)-1)]
      
      
      
      time.serie$Planting   <- start
      time.serie1$Planting  <- start
      time.serie2$Planting  <- start
      
      time.serie$Harvest    <- end
      time.serie1$Harvest   <- end
      time.serie2$Harvest   <- end
      
      
      time.serie$Duration   <- ifelse(test = time.serie$Planting  < time.serie$Harvest, yes = "One year", no = "Two years")
      time.serie1$Duration  <- ifelse(test = time.serie1$Planting < time.serie1$Harvest, yes = "One year", no = "Two years")
      time.serie2$Duration  <- ifelse(test = time.serie2$Planting < time.serie2$Harvest, yes = "One year", no = "Two years")
      
      duration <- time.serie$Duration
      
      time.serie  <- time.serie [, 1:(ncol(time.serie)-3)]
      time.serie1 <- time.serie1[, 1:(ncol(time.serie1)-3)]
      time.serie2 <- time.serie2[, 1:(ncol(time.serie2)-3)]
      
      
      
      
      if(duration== "One year"){
        suppressMessages(library(tidyverse))
        suppressMessages(library(compiler))
        
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        Y <- time.serie1
        Y <- Y %>% gather(key = Date, value = Value, -(cellID:lat))
        Y$Year <- lubridate::year(as.Date(Y$Date))
        Y$Yday <- lubridate::yday(as.Date(Y$Date))
        Y <- Y %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        Z <- time.serie2
        Z <- Z %>% gather(key = Date, value = Value, -(cellID:lat))
        Z$Year <- lubridate::year(as.Date(Z$Date))
        Z$Yday <- lubridate::yday(as.Date(Z$Date))
        Z <- Z %>% group_by(Year) %>% dplyr::filter(Yday >= start& Yday <= end)
        
        
        #Vapour pressure deficit
        calc_vpd <- function(srad, tmin, tmax){
          
          #constants
          albedo <- 0.2
          vpd_cte <- 0.7
          
          #soil heat flux parameters
          a_eslope=611.2
          b_eslope=17.67
          c_eslope=243.5
          
          #input parameters
          tmean <- (tmin+tmax)/2
          
          #net radiation
          rn = (1-albedo) * srad
          
          #soil heat flux
          eslope=a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
          
          #estimate vpd
          esat_min=0.61120*exp((17.67*tmin)/(tmin+243.5))
          esat_max=0.61120*exp((17.67*tmax)/(tmax+243.5))
          vpd=vpd_cte*(esat_max-esat_min) #kPa
          return(vpd)
        }
        
        vpd <- as.data.frame(calc_vpd(srad = Z$Value, tmin = Y$Value, tmax = X$Value ))
        vpd <- data.frame(cellID = X$cellID, lon= X$lon, lat= X$lat, Value= vpd, Year= X$Year, Yday=X$Yday) 
        names(vpd)<- c("cellID","lon","lat","Value","Year","Yday")
        results <- vpd %>% group_by(Year) %>% summarize(median = median(Value)) 
        results <- data.frame(cellID = X$cellID, lon= X$lon, lat= X$lat,results)
        
      }else {
        
        suppressMessages(library(tidyverse))
        suppressMessages(library(compiler))
        ###  
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% filter(Yday %in% c(start:365, 1:end))
        X <- X[-(1:(end)),]
        X <- X[-((nrow(X)-(365-start)): nrow(X)),]
        X$condition <- NA
        for(j in 1:nrow(X)){
          X$condition[j] <- X$Yday[j+1] - X$Yday[j]
        }
        
        c <- unique(X$condition)
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
        X <- na.omit(X)
        #
        
        Y <- time.serie1
        Y <- Y %>% gather(key = Date, value = Value, -(cellID:lat))
        Y$Year <- lubridate::year(as.Date(Y$Date))
        Y$Yday <- lubridate::yday(as.Date(Y$Date))
        Y <- Y %>% filter(Yday %in% c(start:365, 1:end))
        Y <- Y[-(1:(end)),]
        Y <- Y[-((nrow(Y)-(365-start)): nrow(Y)),]
        
        Y$condition <- NA
        for(j in 1:nrow(Y)){
          Y$condition[j] <- Y$Yday[j+1] - Y$Yday[j]
        }
        c<- unique(Y$condition)
        Y$condition[which(is.na(Y$condition))] <- c[3]
        
        chngs <- which(Y$condition == c[3])
        for(j in 1:length(chngs)){
          if(j == 1){
            Y$condition2[1:chngs[j]] <- j
          } else {
            Y$condition2[(chngs[j-1]+1):(chngs[j])] <- j
          }
        }
        Y$condition2[is.na(Y$condition2)] <- length(chngs)
        Y$condition <- NULL
        names(Y)[ncol(Y)] <- "condition"
        rownames(Y)<- 1:nrow(Y)
        
        Y <- na.omit(Y)
        #
        Z <- time.serie2
        Z <- Z %>% gather(key = Date, value = Value, -(cellID:lat))
        Z$Year <- lubridate::year(as.Date(Z$Date))
        Z$Yday <- lubridate::yday(as.Date(Z$Date))
        Z <- Z %>% filter(Yday %in% c(start:365, 1:end))
        Z <- Z[-(1:(end)),]
        Z <- Z[-((nrow(Z)-(365-start)): nrow(Z)),]
        
        Z$condition <- NA
        for(j in 1:nrow(Z)){
          Z$condition[j] <- Z$Yday[j+1] - Z$Yday[j]
        }
        c<- unique(Z$condition)
        Z$condition[which(is.na(Z$condition))] <- c[3]
        
        chngs <- which(Z$condition == c[3])
        for(j in 1:length(chngs)){
          if(j == 1){
            Z$condition2<-NA
            Z$condition2[1:chngs[j]] <- j
          } else {
            Z$condition2[(chngs[j-1]+1):(chngs[j])] <- j
          }
        }
        Z$condition2[is.na(Z$condition2)] <- length(chngs)
        Z$condition <- NULL
        names(Z)[ncol(Z)] <- "condition"
        rownames(Z)<- 1:nrow(Z)
        
        Z <- na.omit(Z) # Dias con tmax > 40 
        
        #Vapour pressure deficit
        calc_vpd <- function(srad, tmin, tmax){
          
          #constants
          albedo <- 0.2
          vpd_cte <- 0.7
          
          #soil heat flux parameters
          a_eslope=611.2
          b_eslope=17.67
          c_eslope=243.5
          
          #input parameters
          tmean <- (tmin+tmax)/2
          
          #net radiation
          rn = (1-albedo) * srad
          
          #soil heat flux
          eslope=a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
          
          #estimate vpd
          esat_min=0.61120*exp((17.67*tmin)/(tmin+243.5))
          esat_max=0.61120*exp((17.67*tmax)/(tmax+243.5))
          vpd=vpd_cte*(esat_max-esat_min) #kPa
          return(vpd)
        }
        
        vpd <- as.data.frame(calc_vpd(srad = Z$Value, tmin = Y$Value, tmax = X$Value ))
        vpd <- data.frame(cellID = X$cellID, lon= X$lon, lat= X$lat, Value= vpd, Year= X$Year, Yday=X$Yday) 
        names(vpd)<- c("cellID","lon","lat","Value","Year","Yday")
        results <- vpd %>% group_by(Year) %>% summarize(median = median(Value)) 
        results <- data.frame(cellID = unique(X$cellID), lon= unique(X$lon), lat=unique(X$lat),results)
        
      }
      return(results)
      
    }, mc.cores = ncores, mc.preschedule = F))
    
    tabla <- do.call(rbind, indexes_been)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    
  }else {
    cat(">>> VPD have been already calculated ...\n")
  }
  
}