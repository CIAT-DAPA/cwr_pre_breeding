# CWR pre-breeding characterising testing environments: calculating general indices
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

generalIndices <- function(crop = "Bean", continent = "Europa"){
  
  output <- paste0(root, "/CWR_pre-breeding/Results/", crop, "/General_indices/", tolower(crop), "_general_indices_", tolower(continent), ".rds")
  if(!file.exists(output)){
    
    # Load climate data
    cat(">>> Starting process for", crop, "in", continent, "continent\n\n")
    cat(">>> Loading climate data ...\n")
    prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_', tolower(continent), '.rds'))
    tmax <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmax/tmax_filtered_', tolower(continent), '.rds'))
    tmin <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_', tolower(continent), '.rds'))
    
    prec$bean_coordinates <- NULL
    tmax$bean_coordinates <- NULL
    tmin$bean_coordinates <- NULL
    
    # Load crop cycle
    cat(">>> Loading crop cycle data ...\n")
    if(crop %in% c("Bean")){
      
      # Planting dates
      planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
      planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
      # Harversting dates
      harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
      harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
      
    } else {
      if(crop == "Sunflower"){
        
        # Planting dates
        planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Sunflower_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
        planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
        # Harversting dates
        harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Sunflower_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
        harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
        
      } else {
        if(crop == "Potato"){
          
          # Planting dates
          planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
          planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
          # Harversting dates
          harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
          harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
          
        } else {
          if(crop == "Finger_millet"){
            
            # Planting dates
            planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Millet_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
            planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
            # Harversting dates
            harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Millet_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
            harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
            
          }
        }
      }
    }
    
    # Extract important dates
    prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")])
    prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")])
    prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")
    
    tmax$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmax[,c("lon", "lat")])
    tmax$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmax[,c("lon", "lat")])
    tmax$Duration <- ifelse(test = tmax$Planting < tmax$Harvest, yes = "One year", no = "Two years")
    
    tmin$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmin[,c("lon", "lat")])
    tmin$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmin[,c("lon", "lat")])
    tmin$Duration <- ifelse(test = tmin$Planting < tmin$Harvest, yes = "One year", no = "Two years")
    rm(planting_rf_ggcmi, harvest_rf_ggcmi)
    
    # Restricting study area to crop area
    cat(">>> Restricting study area to crop area ...\n")
    crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/", crop, "/database/area_base.rds"))
    prec <- dplyr::filter(prec, prec$cellID %in% crop_area$cellID)
    prec <- prec[!is.na(prec$cellID),]
    tmax <- dplyr::filter(tmax, tmax$cellID %in% crop_area$cellID)
    tmax <- tmax[!is.na(tmax$cellID),]
    tmin <- dplyr::filter(tmin, tmin$cellID %in% crop_area$cellID)
    tmin <- tmin[!is.na(tmin$cellID),]
    
    # Define concordant pixels to calculate indices
    pixelList <- Reduce(intersect, list(prec[,'cellID'], tmax[,'cellID'], tmin[,'cellID']))
    
    # Index calculation
    cat(">>> Calculating agroclimatic indices ...\n")
    general_indices <- parallel::mclapply(X = 1:length(pixelList), FUN = function(i){
      
      cat(paste0("Processing pixel: ", pixelList[i], "\n"))
      
      # Define parameters
      duration <- prec$Duration[prec$cellID == pixelList[i]]
      start <- prec$Planting[prec$cellID == pixelList[i]]
      end <- prec$Harvest[prec$cellID == pixelList[i]]
      
      # Calculate indices
      if(duration == "One year"){
        
        # 1. Tmean: Mean daily temperature averaged for a specified period
        time.serie <- (tmax[which(tmax$cellID == pixelList[i]), 1:(ncol(tmax)-3)] + tmin[which(tmin$cellID == pixelList[i]), 1:(ncol(tmin)-3)])/2
        X <- time.serie; rm(time.serie)
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        tmean <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TMEAN = mean(Value, na.rm = T))
        tmean <- tmean %>% as.data.frame
        names(tmean)[2] <- "Value"; tmean$Variable <- "TMEAN"
        
        # 2. GDD_1: Crop duration. Growing degree days calculated using a capped-top function with TB=10 ?C
        calc_cdur <- function(TMEAN, season_ini=1, season_end=365, t_thresh=35){
          tmean <- mean(TMEAN[season_ini:season_end], na.rm=T)
          if (tmean > t_thresh) {cdur <- tmean - t_thresh} else {cdur <- 0}
          return(cdur)
        }
        library(compiler)
        calc_cdurCMP <- cmpfun(calc_cdur)
        gdd1 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(GDD_1 = calc_cdurCMP(Value, t_thresh = 10))
        gdd1 <- gdd1 %>% as.data.frame
        names(gdd1)[2] <- "Value"; gdd1$Variable <- "GDD_1"
        
        # 3. GDD_2: Crop duration. Growing degree days calculated using a capped-top function with TB=25 ?C
        gdd2 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(GDD_2 = calc_cdurCMP(Value, t_thresh = 25))
        gdd2 <- gdd2 %>% as.data.frame
        names(gdd2)[2] <- "Value"; gdd2$Variable <- "GDD_2"
        
        # 4. ND_t35: Heat stress. Total number of days with maximum temperature greater or equal to 35 ?C
        time.serie <- tmax[which(tmax$cellID == pixelList[i]), 1:(ncol(tmax)-3)]
        X <- time.serie; rm(time.serie)
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        calc_hts <- function(TMAX, season_ini=1, season_end=365, t_thresh=35) {
          hts <- length(which(TMAX[season_ini:season_end] >= t_thresh))
          return(hts)
        }
        library(compiler)
        calc_htsCMP <- cmpfun(calc_hts)
        ndt35 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(ND_t35 = calc_htsCMP(Value, t_thresh = 35))
        ndt35 <- ndt35 %>% as.data.frame
        names(ndt35)[2] <- "Value"; ndt35$Variable <- "ND_t35"
        
        # 5. TOTRAIN: Total precipitation
        time.serie <- prec[which(prec$cellID == pixelList[i]), 1:(ncol(prec)-3)]
        X <- time.serie; rm(time.serie)
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
        totrain <- totrain %>% as.data.frame
        names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
        
        # 6. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
        dr_stress <- function(PREC, p_thresh = 1){
          runs <- rle(PREC < p_thresh)
          cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
          return(cons_days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
        cdd <- cdd %>% as.data.frame
        names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
        
        # 7. P5D: Maximum 5-day running average precipitation (Flooding)
        run_avg <- function(x){
          z <- caTools::runmean(x, k = 5, endrule = 'NA')
          z <- max(z, na.rm = TRUE)
          return(z)
        }
        p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
        p5d <- p5d %>% as.data.frame
        names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
        
        # 8. P_95: 95th percentile of daily precipitation (Erosion risk)
        p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
        p_95 <- p_95 %>% as.data.frame
        names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
        
        # Putting all results together
        results <- data.frame(cellID = unique(X$cellID), rbind(tmean, gdd1, gdd2, ndt35, totrain, cdd, p5d, p_95))
        
      } else {
        if(duration == "Two years"){
          
          # 1. Tmean: Mean daily temperature averaged for a specified period
          time.serie <- (tmax[which(tmax$cellID == pixelList[i]), 1:(ncol(tmax)-3)] + tmin[which(tmin$cellID == pixelList[i]), 1:(ncol(tmin)-3)])/2
          X <- time.serie; rm(time.serie)
          X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
          X$Year <- lubridate::year(as.Date(X$Date))
          X$Yday <- lubridate::yday(as.Date(X$Date))
          X <- X %>% filter(Yday %in% c(start:365, 1:end))
          X <- X[-(1:(end)),]
          X <- X[-((nrow(X)-(365-start)): nrow(X)),]
          
          tmean <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TMEAN = mean(Value, na.rm = T))
          tmean <- tmean %>% as.data.frame
          names(tmean)[2] <- "Value"; tmean$Variable <- "TMEAN"
          
          # 2. GDD_1: Crop duration. Growing degree days calculated using a capped-top function with TB=10 ?C
          calc_cdur <- function(TMEAN, season_ini=1, season_end=365, t_thresh=35){
            tmean <- mean(TMEAN[season_ini:season_end], na.rm=T)
            if (tmean > t_thresh) {cdur <- tmean - t_thresh} else {cdur <- 0}
            return(cdur)
          }
          library(compiler)
          calc_cdurCMP <- cmpfun(calc_cdur)
          gdd1 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(GDD_1 = calc_cdurCMP(Value, t_thresh = 10))
          gdd1 <- gdd1 %>% as.data.frame
          names(gdd1)[2] <- "Value"; gdd1$Variable <- "GDD_1"
          
          # 3. GDD_2: Crop duration. Growing degree days calculated using a capped-top function with TB=25 ?C
          gdd2 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(GDD_2 = calc_cdurCMP(Value, t_thresh = 25))
          gdd2 <- gdd2 %>% as.data.frame
          names(gdd2)[2] <- "Value"; gdd2$Variable <- "GDD_2"
          
          # 4. ND_t35: Heat stress. Total number of days with maximum temperature greater or equal to 35 ?C
          time.serie <- tmax[which(tmax$cellID == pixelList[i]), 1:(ncol(tmax)-3)]
          X <- time.serie; rm(time.serie)
          X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
          X$Year <- lubridate::year(as.Date(X$Date))
          X$Yday <- lubridate::yday(as.Date(X$Date))
          X <- X %>% filter(Yday %in% c(start:365, 1:end))
          X <- X[-(1:(end)),]
          X <- X[-((nrow(X)-(365-start)): nrow(X)),]
          
          calc_hts <- function(TMAX, season_ini=1, season_end=365, t_thresh=35) {
            hts <- length(which(TMAX[season_ini:season_end] >= t_thresh))
            return(hts)
          }
          library(compiler)
          calc_htsCMP <- cmpfun(calc_hts)
          ndt35 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(ND_t35 = calc_htsCMP(Value, t_thresh = 35))
          ndt35 <- ndt35 %>% as.data.frame
          names(ndt35)[2] <- "Value"; ndt35$Variable <- "ND_t35"
          
          # 5. TOTRAIN: Total precipitation
          time.serie <- prec[which(prec$cellID == pixelList[i]), 1:(ncol(prec)-3)]
          X <- time.serie; rm(time.serie)
          X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
          X$Year <- lubridate::year(as.Date(X$Date))
          X$Yday <- lubridate::yday(as.Date(X$Date))
          X <- X %>% filter(Yday %in% c(start:365, 1:end))
          X <- X[-(1:(end)),]
          X <- X[-((nrow(X)-(365-start)): nrow(X)),]
          
          totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
          totrain <- totrain %>% as.data.frame
          names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
          
          # 6. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
          dr_stress <- function(PREC, p_thresh = 1){
            runs <- rle(PREC < p_thresh)
            cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
            return(cons_days)
          }
          dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
          cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
          cdd <- cdd %>% as.data.frame
          names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
          
          # 7. P5D: Maximum 5-day running average precipitation (Flooding)
          run_avg <- function(x){
            z <- caTools::runmean(x, k = 5, endrule = 'NA')
            z <- max(z, na.rm = TRUE)
            return(z)
          }
          p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
          p5d <- p5d %>% as.data.frame
          names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
          
          # 8. P_95: 95th percentile of daily precipitation (Erosion risk)
          p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
          p_95 <- p_95 %>% as.data.frame
          names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
          
          # Putting all results together
          results <- data.frame(cellID = unique(X$cellID), rbind(tmean, gdd1, gdd2, ndt35, totrain, cdd, p5d, p_95))
          
        }
      }
      
      return(results)
      
    }, mc.cores = 15, mc.preschedule = F)
    tabla <- do.call(rbind, general_indices)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    
    return(cat("Process done\n"))
    
  } else {
    cat(">>> Agroclimatic indices have been already calculated ...\n")
  }
  
}
