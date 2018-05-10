# CWR pre-breeding characterising testing environments: calculating future water balance
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)

# Load packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})

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
src <- paste0(root, "/CWR_pre-breeding/Scripts")

calc_wat_bal <- function(continent = "Europa", gcm = "gcm1", ncores = 10){
  
  source(paste0(src, "/00_calculating_risk_indices.R"))
  # Load climate data
  prec <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_future_climate/rcp85/", gcm, "/chirps/prec_filtered_", tolower(continent), ".rds"))
  tmax <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_future_climate/rcp85/", gcm, "/agmerra_tmax/tmax_filtered_", tolower(continent), ".rds"))
  tmin <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_future_climate/rcp85/", gcm, "/agmerra_tmin/tmin_filtered_", tolower(continent), ".rds"))
  srad <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_future_climate/rcp85/", gcm, "/agmerra_srad/srad_filtered_", tolower(continent), ".rds"))
  
  prec$bean_coordinates <- NULL
  tmax$bean_coordinates <- NULL
  tmin$bean_coordinates <- NULL
  srad$bean_coordinates <- NULL
  
  pixelList <- Reduce(intersect, list(prec[,'cellID'], tmax[,'cellID'], tmin[,'cellID'], srad[,'cellID']))
  
  cat(">>> Loading soil capacity ...\n")
  if(!file.exists(paste0(root, "/CWR_pre-breeding/Input_data/_soils/Soil_capacity/soil_capacity_", tolower(continent), ".rds"))){
    
    # Load soil data
    l <- raster::stack(list.files(path = paste0(root, "/CWR_pre-breeding/Input_data/_soils/AWCh2"), pattern = "*.tif$", full.names = T))
    r <- resample(l, base, method = "ngb")
    rm(l)
    
    soil_data <- data.frame(raster::xyFromCell(object = base, cell = pixelList), raster::extract(x = r, y = raster::xyFromCell(object = base, cell = pixelList))); rm(r)
    names(soil_data)[1:2] <- c("lon", "lat")
    soil_data$root_depth <- 100
    soil_data <- soil_data %>% select(lon, lat, root_depth, AWCh2_M_sl1_10km_ll:AWCh2_M_sl6_10km_ll)
    
    depths <- c(25, 100, 225, 450, 800, 1500)
    colnames(soil_data)[4:ncol(soil_data)] <- paste0("d.", depths); rm(depths)
    
    # Calculate soil capacity
    soilcap_calc_mod <- function(x, minval, maxval) {
      if(!is.na(x[3])){
        rdepth <- max(c(x[3], minval)) #cross check
        rdepth <- min(c(rdepth, maxval)) #cross-check
        wc_df <- data.frame(depth = c(2.5, 10, 22.5, 45, 80, 150), wc = (x[4:9]) * .01)
        if(!rdepth %in% wc_df$depth) {
          wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
          wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
          y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
          x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
          ya <- (rdepth-x1)/(x2-x1) * (y2-y1) + y1
          wc_df <- rbind(wc_df1, data.frame(depth = rdepth, wc = ya), wc_df2)
        }
        wc_df <- wc_df[which(wc_df$depth <= rdepth),]
        wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
        wc_df$soilcap <- wc_df$soilthick * wc_df$wc
        soilcp <- sum(wc_df$soilcap) * 10 # in mm
        return(soilcp)
      } else {
        soilcp <- NA
        return(soilcp)
      }
    }
    soil_data$soilcp <- apply(soil_data, 1, FUN = soilcap_calc_mod, minval = 45, maxval = 100)
    soil_data$cellID <- cellFromXY(object = base, xy = soil_data[, c("lon", "lat")])
    soil_data <- soil_data %>% select(cellID, lon, lat, root_depth, d.25:soilcp)
    saveRDS(object = soil_data, file = paste0(root, "/CWR_pre-breeding/Input_data/_soils/Soil_capacity/soil_capacity_", tolower(continent), ".rds"))
    
  } else {
    
    if(file.exists(paste0(root, "/CWR_pre-breeding/Input_data/_soils/Soil_capacity/soil_capacity_", tolower(continent), ".rds"))){
      soil_data <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_soils/Soil_capacity/soil_capacity_", tolower(continent), ".rds"))
    } else {
      cat(">>> Please check what happened ...\n")
    }
    
  }
  
  # Calculate water balance per pixel
  daysList <- Reduce(intersect, list(colnames(tmax[,-c(1:3)]), colnames(tmin[,-c(1:3)]),
                                     colnames(prec[,-c(1:3)]), colnames(srad[,-c(1:3)])))
  
  if(OSys == "Linux"){
    
    watbalList <- parallel::mclapply(1:length(pixelList), FUN = function(j){
      
      if(!file.exists(paste0(root,
                             "/CWR_pre-breeding/Input_data/_soils/Water_balance_future/rcp85/",
                             gcm,
                             "/",
                             continent,
                             "/cellID_",
                             pixelList[j],
                             ".rds"))){
        
        out_all <- soil_data[which(soil_data$cellID == pixelList[j]), c('cellID', 'lon', 'lat')]
        out_all <- do.call("rbind", replicate(length(daysList), out_all, simplify = FALSE))
        out_all$SRAD <- as.numeric(srad[which(srad$cellID == pixelList[j]), match(daysList, colnames(srad))])
        out_all$TMIN <- as.numeric(tmin[which(tmin$cellID == pixelList[j]), match(daysList, colnames(tmin))])
        out_all$TMAX <- as.numeric(tmax[which(tmax$cellID == pixelList[j]), match(daysList, colnames(tmax))])
        out_all$RAIN <- as.numeric(prec[which(prec$cellID == pixelList[j]), match(daysList, colnames(prec))])
        rownames(out_all) <- daysList
        soilcp <- soil_data[which(soil_data$cellID == pixelList[j]), 'soilcp']
        watbal_loc <- watbal_wrapper(out_all = out_all, soilcp = soilcp)
        
        general_watbal <- tibble::tibble(cellID = pixelList[j])
        general_watbal <- cbind(general_watbal, raster::xyFromCell(object = base, cell = pixelList[j])) %>% as.tibble()
        colnames(general_watbal)[2:3] <- c("lon", "lat")
        general_watbal$watbal <- watbal_loc %>% list
        
        saveRDS(object = general_watbal,
                file = paste0(root,
                              "/CWR_pre-breeding/Input_data/_soils/Water_balance_future/rcp85/",
                              gcm,
                              "/",
                              continent,
                              "/cellID_",
                              pixelList[j],
                              ".rds"))
        
      } else {
        cat(paste0("Pixel: ", pixelList[j], " has been already processed.\n"))
      }
      
    }, mc.cores = ncores, mc.preschedule = F)
    
  } else {
    if(OSys == "Windows"){
      
      watbalList <- parallelsugar::mclapply(1:length(pixelList), FUN = function(j){
        
        if(!file.exists(paste0(root,
                               "/CWR_pre-breeding/Input_data/_soils/Water_balance_future/rcp85/",
                               gcm,
                               "/",
                               continent,
                               "/cellID_",
                               pixelList[j],
                               ".rds"))){
          
          out_all <- soil_data[which(soil_data$cellID == pixelList[j]), c('cellID', 'lon', 'lat')]
          out_all <- do.call("rbind", replicate(length(daysList), out_all, simplify = FALSE))
          out_all$SRAD <- as.numeric(srad[which(srad$cellID == pixelList[j]), match(daysList, colnames(srad))])
          out_all$TMIN <- as.numeric(tmin[which(tmin$cellID == pixelList[j]), match(daysList, colnames(tmin))])
          out_all$TMAX <- as.numeric(tmax[which(tmax$cellID == pixelList[j]), match(daysList, colnames(tmax))])
          out_all$RAIN <- as.numeric(prec[which(prec$cellID == pixelList[j]), match(daysList, colnames(prec))])
          rownames(out_all) <- daysList
          soilcp <- soil_data[which(soil_data$cellID == pixelList[j]), 'soilcp']
          watbal_loc <- watbal_wrapper(out_all = out_all, soilcp = soilcp)
          
          general_watbal <- tibble::tibble(cellID = pixelList[j])
          general_watbal <- cbind(general_watbal, raster::xyFromCell(object = base, cell = pixelList[j])) %>% as.tibble()
          colnames(general_watbal)[2:3] <- c("lon", "lat")
          general_watbal$watbal <- watbal_loc %>% list
          
          saveRDS(object = general_watbal,
                  file = paste0(root,
                                "/CWR_pre-breeding/Input_data/_soils/Water_balance_future/rcp85/",
                                gcm,
                                "/",
                                continent,
                                "/cellID_",
                                pixelList[j],
                                ".rds"))
          
        } else {
          cat(paste0("Pixel: ", pixelList[j], " has been already processed.\n"))
        }
        
      }, mc.cores = ncores, mc.preschedule = F)
      
    }
  }
  
  return(cat(">>> Process done\n"))
  
}

calc_wat_bal(continent = "Europa", gcm = "gcm1", ncores = 7) # Climate
calc_wat_bal(continent = "Oceania", gcm = "gcm1", ncores = 10) # Africa
