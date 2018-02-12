### CWR_Pre_breeding
### General Index 
### load Packages 

suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
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
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
suppressMessages(if(!require(googlesheets)){install.packages("googlesheets");library(googlesheets)}else{library(googlesheets)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})

# Path options
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

prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_africa.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_europa.rds'))
# prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_asia.rds'))
#prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_oceania.rds'))
prec$bean_coordinates <- NULL

# Planting dates
planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(planting_rf_ggcmi)
prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(harvest_rf_ggcmi)
prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")

crop  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Precense_data/area_base.rds"))
prec <- dplyr::filter(prec, prec$cellID %in% crop$cellID)
prec<- prec[!is.na(prec$cellID),]

# Just for verifying
prec[1:5,(ncol(prec)-5):ncol(prec)]
table(prec$Duration)

library(parallel)
system.time( prec_indexes <- mclapply(1:nrow(prec), function(i){
  cat(paste0("Processed pixel:", i, "\n"))
  # Parameters
  duration <- prec$Duration[i]
  start <- prec$Planting[i]
  end <- prec$Harvest[i]
  
  # Just one pixel
  time.serie <- prec[i, 1:(ncol(prec)-3)]
  
  if(duration == "One year"){
    
    suppressMessages(library(tidyverse))
    suppressMessages(library(compiler))
    
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
    
    # Total precipitation
    totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
    totrain <- totrain %>% as.data.frame
    names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
    
    # Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
    dr_stress <- function(PREC, p_thresh = 1){
      runs <- rle(PREC < p_thresh)
      cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
    cdd <- cdd %>% as.data.frame
    names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
    
    # # Flooding: Maximum 5-day running average precipitation
    require(caTools)
    run_avg <- function(x){
      z <- caTools::runmean(x, k = 5, endrule = 'NA')
      z <- max(z, na.rm = TRUE)
      return(z)
    }
    p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
    p5d <- p5d %>% as.data.frame
    names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
    
    
    # Erosion risk: 95th percentile of daily precipitation
    p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
    p_95 <- p_95 %>% as.data.frame
    names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
    
    results <- data.frame(cellID = unique(X$cellID), rbind(totrain, cdd,p5d, p_95))
    
  } else {
    suppressMessages(library(tidyverse))
    suppressMessages(library(compiler))
    X <- time.serie
    X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
    X$Year <- lubridate::year(as.Date(X$Date))
    X$Yday <- lubridate::yday(as.Date(X$Date))
    X <- X %>% filter(Yday %in% c(start:365, 1:end))
    X <- X[-(1:(end)),]
    X <- X[-((nrow(X)-(365-start)): nrow(X)),]
    
    # Total precipitation
    totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
    totrain <- totrain %>% as.data.frame
    names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
    
    # Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
    dr_stress <- function(PREC, p_thresh = 1){
      runs <- rle(PREC < p_thresh)
      cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      return(cons_days)
    }
    dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
    cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
    cdd <- cdd %>% as.data.frame
    names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
    
    require(caTools)
    # Flooding: Maximum 5-day running average precipitation
    run_avg <- function(x){
      z <- caTools::runmean(x, k = 5, endrule = 'NA')
      z <- max(z, na.rm = TRUE)
      return(z)
    }
    p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
    p5d <- p5d %>% as.data.frame
    names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
    
    # Erosion risk: 95th percentile of daily precipitation
    p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
    p_95 <- p_95 %>% as.data.frame
    names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
    
    results <- data.frame(cellID = unique(X$cellID), rbind(totrain, cdd, p5d, p_95))
  }
  return(results)
}, mc.cores = 15, mc.preschedule = F))
tabla <- do.call(rbind,prec_indexes)

saveRDS(tabla,paste0(root,"/CWR_pre-breeding/Input_data/general_indexes/bean_index_tabla_rbind_america.rds"))
# saveRDS(tabla,paste0(root,"/CWR_pre-breeding/Input_data/general_indexes/bean_index_tabla_rbind_africa.rds"))
# saveRDS(tabla,paste0(root,"/CWR_pre-breeding/Input_data/general_indexes/bean_index_tabla_rbind_europa.rds"))
# saveRDS(tabla,paste0(root,"/CWR_pre-breeding/Input_data/general_indexes/bean_index_tabla_rbind_asia.rds"))
# saveRDS(tabla,paste0(root,"/CWR_pre-breeding/Input_data/general_indexes/bean_index_tabla_rbind_oceania.rds"))




