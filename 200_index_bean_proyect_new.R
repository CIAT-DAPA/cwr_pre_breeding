#CWR pre-breeding characterising testing environments: calculating general indices
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

beanIndices <- function(continent = "Europa"){
  
  output <- paste0(root, "/CWR_pre-breeding/Proyect_bean/Calculate_index_heat/bean_crop_indices_", tolower(continent), ".rds")
  if(!file.exists(output)){
    
    
    # Load climate data
    cat(">>> Starting process for bean in", continent, "continent\n\n")
    cat(">>> Loading climate data ...\n")
    
    
    tmax <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmax/tmax_filtered_', tolower(continent), '.rds'))
    tmin <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_', tolower(continent), '.rds')) 
    
    tmax$bean_coordinates <- NULL
    tmin$bean_coordinates <- NULL
    
    # Load crop cycle
    cat(">>> Loading crop cycle data ...\n")
    # Planting dates
    planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
    planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
    # Harversting dates
    harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
    harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
    
    # Restricting study area to crop area
    cat(">>> Restricting study area to crop area ...\n")
    crop<- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Bean/database/area_base.rds"))
    
    cell <- crop$cellID
    tmax <- filter(tmax,  cellID  %in% cell)
    tmin <- filter (tmin, cellID %in% cell)
    tmean<-((tmax[,4:(ncol(tmax)-1)])+(tmin[,4:(ncol(tmin)-1)]))*0.5
    tmean<- data.frame(cellID= tmax$cellID, lon= tmax$lon, lat= tmax$lat, tmean)
    names(tmean)[4:ncol(tmean)] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(tmean)[4:ncol(tmean)]), fixed = T))
    rm(crop)
    
    tmax  <- tmax[!duplicated(tmax$cellID),]
    tmin  <- tmin[!duplicated(tmin$cellID),]
    tmean <- tmean[!duplicated(tmean$cellID),]
    
    
    # Extract important dates
    tmax$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmax[,c("lon", "lat")]);
    tmax$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmax[,c("lon", "lat")]); 
    tmax$Duration <- ifelse(test = tmax$Planting < tmax$Harvest, yes = "One year", no = "Two years")
    
    tmin$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmin[,c("lon", "lat")]);
    tmin$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmin[,c("lon", "lat")]); 
    tmin$Duration <- ifelse(test = tmin$Planting < tmin$Harvest, yes = "One year", no = "Two years")
    tmin$length <- abs(tmin$Planting - tmin$Harvest)
    
    tmean$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmean[,c("lon", "lat")]);
    tmean$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmean[,c("lon", "lat")]); 
    tmean$Duration <- ifelse(test = tmean$Planting < tmean$Harvest, yes = "One year", no = "Two years")
    
    
    test<- tmin[which(tmin$Duration =="Two years"), ]
    library(doSNOW)
    library(foreach)
    library(parallel)
    library(doParallel)
    
    cores<- detectCores()
    cl<- makeCluster(cores-18)
    registerDoParallel(cl) 
    
    system.time(indexes_been <- foreach(i=1:nrow(tmin)) %dopar% {  ### nrow(tmax)
      tryCatch({ mylist <- list()
      cat(paste0("Processed pixel:", i, "\n"))
      # Parameters
      duration <- tmax$Duration[i]
      start <- tmax$Planting[i]
      end <- tmax$Harvest[i]
      
      # Just one pixel
      time.serie <- tmax[i, 1:(ncol(tmax)-3)]
      time.serie1 <- tmin[i, 1:(ncol(tmin)-4)]
      time.serie2 <- tmean[i, 1:(ncol(tmean)-3)]
      
      if(duration== "One year"){
        suppressMessages(library(tidyverse))
        suppressMessages(library(compiler))
        
        ##NOTA :  X = Tmax , Y= Tmin , Z=Tmean
        
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
        
        # Dias con tmax > 40 
        dr_stress<- function(TMAX, p_thresh = 40){
          days <- sum(TMAX > 40)
          return(days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        Texteme <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(Texteme = dr_stressCMP(Value, p_thresh = 40))
        Texteme <- Texteme %>% as.data.frame
        names(Texteme)[2] <- "Value"; Texteme$Variable <- "Texteme"
        
        # Dias con tmin > 22
        dr_stress<- function(TMIN, p_thresh = 22){
          days <- sum(TMIN > 22)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        days_tmin22 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tmin22 = dr_stressCMP(Value, p_thresh = 22))
        days_tmin22 <- days_tmin22 %>% as.data.frame
        names(days_tmin22)[2] <- "Value"; days_tmin22$Variable <- "days_tmin22"
        
        
        # Dias con tmin > 24
        dr_stress<- function(TMIN, p_thresh = 24){
          days <- sum(TMIN > 24)
          return(days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        days_tmin24 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tmin24 = dr_stressCMP(Value, p_thresh = 24))
        days_tmin24 <- days_tmin24 %>% as.data.frame
        names(days_tmin24)[2] <- "Value"; days_tmin24$Variable <- "days_tmin24"
        
        # Index de calor22 (%)
        long <- tmin$length[i]
        dr_stress<- function(TMIN, p_thresh= 22){
          fic <- (long -30)
          sum_day <- sum(TMIN > p_thresh)
          c <- (sum_day)/fic
          return(c)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        calor22 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor22 = dr_stressCMP(Value, p_thresh = 22))
        calor22 <- calor22 %>% as.data.frame
        names(calor22)[2] <- "Value"; calor22$Variable <- "calor22"
        
        # Index de calor24 (%)
        long <- tmin$length[i]
        dr_stress<- function(TMIN, p_thresh= 22){
          fic <- (long -30)
          sum_day <- sum(TMIN > p_thresh)
          c <- (sum_day)/fic
          return(c)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        calor24 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor24 = dr_stressCMP(Value, p_thresh = 24))
        calor24 <- calor24 %>% as.data.frame
        names(calor24)[2] <- "Value"; calor24$Variable <- "calor24"
        
        # 20<tmean<29
        dr_stress <- function(TMEAN, MIN= 20, MAX= 29 ){
          runs <- rle(TMEAN > MIN & TMEAN <MAX)
          if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
            cons_days <- 0
          } else {
            cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
          }
          return(cons_days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        optimo <- Z %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(optimo = dr_stressCMP(Value, MIN = 20,MAX=27))
        optimo <- optimo %>% as.data.frame
        names(optimo)[2] <- "Value"; optimo$Variable <- "optimo"
        
        
        results<- data.frame(cellID= unique(Y$cellID), rbind(Texteme,days_tmin22,days_tmin24,calor22,calor24,optimo))
        
        mylist[[i]] <- results
        
      }else {
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
        dr_stress<- function(TMAX, p_thresh = 40){
          days <- sum(TMAX > 40)
          return(days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        Texteme <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(Texteme = dr_stressCMP(Value, p_thresh = 40))
        Texteme <- Texteme %>% as.data.frame
        names(Texteme)[2] <- "Value"; Texteme$Variable <- "Texteme"
        
        # Dias con tmin > 22 
        dr_stress<- function(TMIN, p_thresh = 22){
          days <- sum(TMIN > 22)
          return(days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        days_tmin22 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tmin22 = dr_stressCMP(Value, p_thresh = 22))
        days_tmin22 <- days_tmin22 %>% as.data.frame
        names(days_tmin22)[2] <- "Value"; days_tmin22$Variable <- "days_tmin22"
        
        
        # Dias con tmin > 24 
        dr_stress<- function(TMIN, p_thresh = 24){
          days <- sum(TMIN > 24)
          return(days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        days_tmin24 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tmin24 = dr_stressCMP(Value, p_thresh = 24))
        days_tmin24 <- days_tmin24 %>% as.data.frame
        names(days_tmin24)[2] <- "Value"; days_tmin24$Variable <- "days_tmin24"
        
        # # Index de calor22 (%)
        long <- tmin$length[i]
        dr_stress<- function(TMIN, p_thresh= 22){
          fic <- (long -30)
          sum_day <- sum(TMIN > p_thresh)
          c <- (sum_day)/fic
          return(c)
        }
        
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        calor22 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor22 = dr_stressCMP(Value, p_thresh = 22))
        calor22 <- calor22 %>% as.data.frame
        names(calor22)[2] <- "Value"; calor22$Variable <- "calor22"
        # # Index de calor24 (%)
        long <- tmin$length[i]
        dr_stress<- function(TMIN, p_thresh= 24,long1){
          fic <- (long - 30)
          sum_day <- sum(TMIN > p_thresh)
          d <- (sum_day)/fic
          return(d)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        calor24 <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor24 = dr_stressCMP(Value, p_thresh = 24))
        calor24 <- calor24 %>% as.data.frame
        names(calor24)[2] <- "Value"; calor24$Variable <- "calor24"
        
        ## 20<tmean<29
        dr_stress <- function(TMEAN, MIN= 20, MAX= 29 ){
          runs <- rle(TMEAN > MIN & TMEAN <MAX)
          if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
            cons_days <- 0
          } else {
            cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
          }
          return(cons_days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        optimo <- Z %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(optimo = dr_stressCMP(Value, MIN = 20,MAX=27))
        optimo <- optimo %>% as.data.frame
        names(optimo)[2] <- "Value"; optimo$Variable <- "optimo"
        results<- data.frame(cellID= unique(Y$cellID), rbind(Texteme,days_tmin22,days_tmin24,calor22,calor24,optimo))
      }
      mylist[[i]] <- results
      
      }, error=function(e){})
    })
    stopCluster(cl)
    tabla <- do.call(rbind, indexes_been)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    return(cat("Process done\n"))
    
  } else {
    cat(">>> Agroclimatic indices have been already calculated ...\n")
  }
  
}

system.time(beanIndices(continent = "Oceania"))
