# CWR pre-breeding characterising testing environments: calculating similarity index for general indices
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

# Load packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(parallel)){install.packages('parallel'); library(parallel)} else {library(parallel)})
suppressMessages(if(!require(googlesheets)){install.packages("googlesheets");library(googlesheets)}else{library(googlesheets)})
suppressMessages(if(!require(parallelDist)){install.packages("parallelDist");library(parallelDist)}else{library(parallelDist)})

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

calc_similarity <- function(crop = "Bean", ncores = 15){
  
  if(!file.exists()){
    
    # Loading information from Google Drive
    cat(">>> Loading data from Google Drive ...\n")
    trials <- gs_ls("crop_information_states")
    trials <- gs_title("crop_information_states")
    trials %>% gs_browse(ws = "Trials")
    trials <- trials %>% gs_read(ws = "Trials", col_types = cols(Longitude = col_character(),
                                                                 Latitude = col_character()))
    trials$Longitude <- gsub(pattern = ",", replacement = ".", x = trials$Longitude) %>% as.character() %>% as.numeric()
    trials$Latitude <- gsub(pattern = ",", replacement = ".", x = trials$Latitude) %>% as.character() %>% as.numeric()
    
    trials <- trials %>% filter(Crop == crop)
    
    # Loading agroclimatic indices
    ind_africa  <- readRDS(paste0(root, "..."))
    ind_america <- readRDS(paste0(root, "..."))
    ind_asia    <- readRDS(paste0(root, "..."))
    ind_europa  <- readRDS(paste0(root, "..."))
    ind_oceania <- readRDS(paste0(root, "..."))
    
    # Putting all together
    giWorld <- rbind(ind_africa, ind_america, ind_asia, ind_europa, ind_oceania)
    rm(ind_africa, ind_america, ind_asia, ind_europa, ind_oceania)
    
    for(i in 1:nrow(trials)){
      
      # Varying by coordinate
      refCoord <- raster::cellFromXY(object = base, xy = trials[i, c("Longitude", "Latitude")])
      
      # List of pixels
      pixelList <- giWorld$cellID %>% unique %>% sort
      
      # Identify indices from reference pixel
      refIndices <- giWorld %>% filter(cellID == refCoord)
      refIndices <- refIndices %>% spread(key = Year, value = Value)
      
      for(j in 1:length(pixelList)){
        
        # Identify indices from compared pixel
        cmpIndices <- giWorld %>% filter(cellID == pixelList[j])
        cmpIndices <- cmpIndices %>% spread(key = Year, value = Value)
        
        colnames(refIndices)[3:ncol(refIndices)] <- paste0("Y", colnames(refIndices)[3:ncol(refIndices)])
        colnames(cmpIndices)[3:ncol(cmpIndices)] <- paste0("Y", colnames(cmpIndices)[3:ncol(cmpIndices)])
        
        matrix.list <- list(data.matrix(refIndices[,3:ncol(refIndices)]), data.matrix(cmpIndices[,3:ncol(cmpIndices)]))
        
        parallelDist::parDist(x = matrix.list, method = "dtw")
        
      }
      
    }
    
    return()
  } else {
    
  }
  
}