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

## INDEX

#### Data cicle from crop
# Planting dates
planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]


### Load date
tmax <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_america.rds"))
tmin <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_america.rds"))
prec <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds"))








fil <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Bean/database/occ_data_complete_case.csv")
fil <- data.frame(cellId = fil$cellID, lon = fil$lon , lat = fil$lat)
cell <- c (cellFromXY(base, cbind(x=fil$lon, y= fil$lat)))
tmax <- filter(tmax,  cellID  %in% cell)
tmin <- filter (tmin, cellID %in% cell)
prec <- filter (prec, cellID %in% cell)
tmean<-((tmax[,4:(ncol(tmax)-3)])+(tmin[,4:(ncol(tmin)-3)]))*0.5
tmean<- data.frame(cellID= tmax$cellID, lon= tmax$lon, lat= tmax$lat, tmean)
names(tmean)[4:ncol(tmean)] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(tmean)[4:ncol(tmean)]), fixed = T))

rm(fil)



tmax$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmax[,c("lon", "lat")]);
tmax$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmax[,c("lon", "lat")]); 
tmax$Duration <- ifelse(test = tmax$Planting < tmax$Harvest, yes = "One year", no = "Two years")

tmin$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmin[,c("lon", "lat")]);
tmin$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmin[,c("lon", "lat")]); 
tmin$Duration <- ifelse(test = tmin$Planting < tmin$Harvest, yes = "One year", no = "Two years")

tmean$Planting <- raster::extract(x = planting_rf_ggcmi, y = tmean[,c("lon", "lat")]);
tmean$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = tmean[,c("lon", "lat")]); 
tmean$Duration <- ifelse(test = tmean$Planting < tmean$Harvest, yes = "One year", no = "Two years")


# Just for verifying
tmax[1:5,(ncol(tmax)-5):ncol(tmax)]
table(tmax$Duration)

tmin[1:5,(ncol(tmin)-5):ncol(tmin)]
table(tmin$Duration)

tmean[1:5,(ncol(tmean)-5):ncol(tmean)]
table(tmean$Duration)

library(parallel)
system.time( tmax_indexes <- mclapply(1:nrow(tmax), function(i){
  cat(paste0("Processed pixel:", i, "\n"))
  # Parameters
  duration <- tmax$Duration[i]
  start <- tmax$Planting[i]
  end <- tmax$Harvest[i]
  
  # Just one pixel
  time.serie <- tmax[i, 1:(ncol(tmax)-3)]
  time.serie1 <- tmin[i, 1:(ncol(tmin)-3)]
  time.serie2 <- tmean[i, 1:(ncol(tmean)-3)]
  
  if(duration == "One year"){
    suppressMessages(library(tidyverse))
    suppressMessages(library(compiler))
    
    ##NOTA :  X = Tmax , Y= Tmin , Z=Tmean
    
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= start + 90)
    
    
    Y <- time.serie1
    Y <- Y %>% gather(key = Date, value = Value, -(cellID:lat))
    Y$Year <- lubridate::year(as.Date(Y$Date))
    Y$Yday <- lubridate::yday(as.Date(Y$Date))
    Y <- Y %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= start + 90)
    
    Z <- time.serie2
    Z <- Z %>% gather(key = Date, value = Value, -(cellID:lat))
    Z$Year <- lubridate::year(as.Date(Z$Date))
    Z$Yday <- lubridate::yday(as.Date(Z$Date))
    Z <- Z %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= start + 90)
    
    
    # Dias consecutivos con tmax > 40  Texteme
    dr_stress <- function(TMAX, p_thresh = 40){
      runs <- rle(TMAX > p_thresh)
      if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
        cons_days <- 0
      } else {
        cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      }
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    Texteme <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(Texteme = dr_stressCMP(Value, p_thresh = 40))
    Texteme <- Texteme %>% as.data.frame
    names(Texteme)[2] <- "Value"; Texteme$Variable <- "Texteme"
    
    
    # Dias consecutivos con tmin > 27
    dr_stress <- function(TMIN, p_thresh = 27){
      runs <- rle(TMIN > p_thresh)
      if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
        cons_days <- 0
      } else {
        cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      }
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    tmin_ex <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(tmin_ex = dr_stressCMP(Value, p_thresh = 27))
    tmin_ex <- tmin_ex %>% as.data.frame
    names(tmin_ex)[2] <- "Value"; tmin_ex$Variable <- "tmin_ex"
    
    
    # Dias con tmin > 27 
    dr_stress<- function(TMIN, p_thresh = 27){
      days <- sum(TMIN > 27)
      return(days)
    }
    
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    days_tex <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tex = dr_stressCMP(Value, p_thresh = 27))
    days_tex <- days_tex %>% as.data.frame
    names(days_tex)[2] <- "Value"; days_tex$Variable <- "days_tex"
    
    
    # Index de calor (%)
    
    dr_stress<- function(TMIN, p_thresh= 27){
      fic <- 91
      sum_day <- sum(TMIN > p_thresh)
      c <- (fic - sum_day)/fic
      return(c)  
    }
    
    
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    calor <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor = dr_stressCMP(Value, p_thresh = 27))
    calor <- calor %>% as.data.frame
    names(calor)[2] <- "Value"; calor$Variable <- "calor"
    
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
    names(optimo)[2] <- "Value"; optimo$Variable <- "calor"
    
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
    chngs <- which(X$condition == "213")
    for(j in 1:length(chngs)){
      if(j == 1){
        X$condition2[1:chngs[j]] <- j
      } else {
        X$condition2[(chngs[j-1]+1):(chngs[j])] <- j
      }
    }
    X$condition2[is.na(X$condition2)] <- length(chngs)+1
    X$condition <- NULL
    names(X)[ncol(X)] <- "condition"
    X$condition[4285:nrow(X)][] <- 29
    
    corte <- c(1:29)
    lista <- lapply(1:length(corte), function(i){
      tablas <- filter (X, condition== corte[i]) 
      return(tablas)  
    })
    
    tab <- lapply(1:length(lista), function(i){
      tabla <-lista[[i]][1:91,]
      return(tabla)
    })
    X <- do.call(rbind,tab)
    #####
    
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
    chngs <- which(Y$condition == "213")
    for(j in 1:length(chngs)){
      if(j == 1){
        Y$condition2[1:chngs[j]] <- j
      } else {
        Y$condition2[(chngs[j-1]+1):(chngs[j])] <- j
      }
    }
    Y$condition2[is.na(X$condition2)] <- length(chngs)+1
    Y$condition <- NULL
    names(Y)[ncol(Y)] <- "condition"
    Y$condition[4285:nrow(Y)][] <- 29
    
    corte <- c(1:29)
    lista <- lapply(1:length(corte), function(i){
      tablas <- filter (Y, condition== corte[i]) 
      return(tablas)  
    })
    
    tab <- lapply(1:length(lista), function(i){
      tabla <-lista[[i]][1:91,]
      return(tabla)
    })
    Y <- do.call(rbind,tab)
    #####
    
    #####
    
    Z <- time.serie1
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
    chngs <- which(Z$condition == "213")
    for(j in 1:length(chngs)){
      if(j == 1){
        Z$condition2[1:chngs[j]] <- j
      } else {
        Z$condition2[(chngs[j-1]+1):(chngs[j])] <- j
      }
    }
    Z$condition2[is.na(Z$condition2)] <- length(chngs)+1
    Z$condition <- NULL
    names(Z)[ncol(Z)] <- "condition"
    Z$condition[4285:nrow(Z)][] <- 29
    
    corte <- c(1:29)
    lista <- lapply(1:length(corte), function(i){
      tablas <- filter (Z, condition== corte[i]) 
      return(tablas)  
    })
    
    tab <- lapply(1:length(lista), function(i){
      tabla <-lista[[i]][1:91,]
      return(tabla)
    })
    Z <- do.call(rbind,tab)
    
    
    # Dias consecutivos con tmax > 40  Texteme
    dr_stress <- function(TMAX, p_thresh = 40){
      runs <- rle(TMAX > p_thresh)
      if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
        cons_days <- 0
      } else {
        cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      }
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    Texteme <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(Texteme = dr_stressCMP(Value, p_thresh = 40))
    Texteme <- Texteme %>% as.data.frame
    names(Texteme)[2] <- "Value"; Texteme$Variable <- "Texteme"
    
    
    # Dias consecutivos con tmin > 27
    dr_stress <- function(TMIN, p_thresh = 27){
      runs <- rle(TMIN > p_thresh)
      if(max(runs$lengths[runs$values==1], na.rm=TRUE) == -Inf){
        cons_days <- 0
      } else {
        cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      }
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    tmin_ex <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(tmin_ex = dr_stressCMP(Value, p_thresh = 27))
    tmin_ex <- tmin_ex %>% as.data.frame
    names(tmin_ex)[2] <- "Value"; tmin_ex$Variable <- "tmin_ex"
    
    
    # Dias con tmin > 27 
    dr_stress<- function(TMIN, p_thresh = 27){
      days <- sum(TMIN > 27)
      return(days)
    }
    
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    days_tex <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(days_tex = dr_stressCMP(Value, p_thresh = 27))
    days_tex <- days_tex %>% as.data.frame
    names(days_tex)[2] <- "Value"; days_tex$Variable <- "days_tex"
    
    
    # Index de calor (%)
    
    dr_stress<- function(TMIN, p_thresh= 27){
      fic <- 91
      sum_day <- sum(TMIN > p_thresh)
      c <- (fic - sum_day)/fic
      return(c)  
    }
    
    
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    calor <- Y %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(calor = dr_stressCMP(Value, p_thresh = 27))
    calor <- calor %>% as.data.frame
    names(calor)[2] <- "Value"; calor$Variable <- "calor"
    
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
    names(optimo)[2] <- "Value"; optimo$Variable <- "calor"
    
    #####
    
    
    
    
  }
  
}))

