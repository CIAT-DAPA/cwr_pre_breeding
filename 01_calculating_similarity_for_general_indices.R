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

calc_similarity <- function(crop = "Bean", scaled = T, ncores = 15){
  
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
    ind_africa  <- readRDS(paste0(root, "/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_africa.rds"))
    ind_america <- readRDS(paste0(root, "/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_america.rds"))
    ind_asia    <- readRDS(paste0(root, "/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_asia.rds"))
    ind_europa  <- readRDS(paste0(root, "/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_europa.rds"))
    ind_oceania <- readRDS(paste0(root, "/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_oceania.rds"))
    
    # Putting all together
    giWorld <- rbind(ind_africa, ind_america, ind_asia, ind_europa, ind_oceania)
    rm(ind_africa, ind_america, ind_asia, ind_europa, ind_oceania)
    
    for(i in 1:nrow(trials)){
      
      # Varying by coordinate
      refCoord <- raster::cellFromXY(object = base, xy = trials[i, c("Longitude", "Latitude")] %>% as.data.frame)
      
      # List of pixels
      pixelList <- giWorld$cellID %>% unique %>% sort
      
      # Identify indices from reference pixel
      refIndices <- giWorld %>% filter(cellID == refCoord)
      refIndices <- refIndices %>% spread(key = Year, value = Value)
      
      # Calculate climatic similarity between pixels
      similarity <- parallel::mclapply(X = 1:length(pixelList), FUN = function(j){
        
        # Identify indices from compared pixel
        cmpIndices <- giWorld %>% filter(cellID == pixelList[j])
        cmpIndices <- cmpIndices %>% spread(key = Year, value = Value)
        
        colnames(refIndices)[3:ncol(refIndices)] <- paste0("Y", colnames(refIndices)[3:ncol(refIndices)])
        colnames(cmpIndices)[3:ncol(cmpIndices)] <- paste0("Y", colnames(cmpIndices)[3:ncol(cmpIndices)])
        
        matrix.list <- list(data.matrix(refIndices[,3:ncol(refIndices)]), data.matrix(cmpIndices[,3:ncol(cmpIndices)]))
        
        # Calculate similarity through multivariate DTW
        calc_multivariate_dtw <- function(data = matrix.list, scaled = T){
          if(scaled){
            means_list1 <- apply(X = matrix.list[[1]], MARGIN = 1, FUN = mean)
            means_list2 <- apply(X = matrix.list[[2]], MARGIN = 1, FUN = mean)
            sds_list1 <- apply(X = matrix.list[[1]], MARGIN = 1, FUN = sd)
            sds_list2 <- apply(X = matrix.list[[1]], MARGIN = 1, FUN = sd)
            matrix.list[[1]] <- t(apply(X = matrix.list[[1]], MARGIN = 1, FUN = scale))
            matrix.list[[2]] <- t(apply(X = matrix.list[[2]], MARGIN = 1, FUN = scale))
            matrix.list[[1]][which(sds_list1 == 0),] <- means_list1[which(sds_list1 == 0)]
            matrix.list[[2]][which(sds_list2 == 0),] <- means_list2[which(sds_list2 == 0)]
            rm(means_list1, means_list2, sds_list1, sds_list2)
            calcDTW <- parallelDist::parDist(x = matrix.list, method = "dtw")
          } else {
            calcDTW <- parallelDist::parDist(x = matrix.list, method = "dtw")
          }
          return(calcDTW)
        }
        dtwObj <- calc_multivariate_dtw(data = matrix.list, scaled = scaled)
        results <- data.frame(cellID = pixelList[j], dtw = dtwObj)
        return(results)
        
      })
      similarity <- do.call(rbind, similarity)
      
      # Save a table with the following columns: DTW calculated, DTW categorized
      
    }
    
    return()
  } else {
    
  }
  
}