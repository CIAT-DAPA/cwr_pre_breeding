# CWR pre-breeding characterising testing environments: calculating current water balance
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
suppressMessages(if(!require(doSNOW)){install.packages('doSNOW'); library(doSNOW)} else {library(doSNOW)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} else {library(doParallel)})

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

continentList <- c("Oceania") # "Europa", "Oceania", "America", "Asia", "Africa"
lapply(1:length(continentList), function(j){
  
  wb_tbl <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_soils/Water_balance/watbal_", tolower(continentList[j]), ".rds"))
  cellID <- wb_tbl$cellID %>% as.character
  
  library(doSNOW)
  library(foreach)
  library(parallel)
  library(doParallel)
  
  cores <- detectCores()
  
  if(OSys == "Linux"){
    cl <- parallel::makeCluster(cores - 3)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- parallel::makeCluster(cores - 3)
    doSNOW::registerDoSNOW(cl)
  }
  
  clusterExport(cl, c('root', 'continentList')) # , 'wb_tbl', 'cellID'
  tmp <- foreach(i = 1:nrow(wb_tbl)) %dopar% { 
    df <- wb_tbl[i,]
    outfile <- paste0(root, "/CWR_pre-breeding/Input_data/_soils/Water_balance/", continentList[j], "/cellID_", cellID[i], ".rds")
    if(!file.exists(outfile)){
      saveRDS(df, file = outfile)
    }
    
  }
  
  parallel::stopCluster(cl)
  
  rm(wb_tbl); rm(cellID)
  rm(cores); rm(cl); rm(tmp)
  g <- gc(); rm(g)
  
  return(cat("Process done for continent", continentList[j], "!\n"))
  
})
