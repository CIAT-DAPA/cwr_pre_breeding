# CWR pre-breeding characterising testing environments: calculating current water balance
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

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
}; rm(OSys)

calc_wat_bal <- function(continent = "Europa"){
  
  # Load climate data
  prec <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_current_climate/chirps/prec_filtered_", tolower(continent), ".rds"))
  tmax <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmax/tmax_filtered_", tolower(continent), ".rds"))
  tmin <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_", tolower(continent), ".rds"))
  srad <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_current_climate/agmerra_srad/srad_filtered_", tolower(continent), ".rds"))
  
  prec$bean_coordinates <- NULL
  tmax$bean_coordinates <- NULL
  tmin$bean_coordinates <- NULL
  srad$bean_coordinates <- NULL
  
  pixelList <- Reduce(intersect, list(prec[,'cellID'], tmax[,'cellID'], tmin[,'cellID'], srad[,'cellID']))
  
  # Load soil data
  l <- raster::stack(list.files(path = paste0(root, "/CWR_pre-breeding/Input_data/_soils/AWCh2"), pattern = "*.tif$", full.names = T))
  r <- resample(l, base, method = "ngb")
  rm(l)
  # Next steps: Obtain coordinates from pixelList to extract depth values from raster
  
  # Calculate soil capacity
  
  
  # Calculate water balance per pixel
  
  
  return(cat(">>> Process done\n"))
  
}