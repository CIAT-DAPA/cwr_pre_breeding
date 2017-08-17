# CWR pre-breeding characterising testing environments: Identify crop cycle per pixel
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#       Load packages           #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
# suppressMessages(if(!require(doMC)){install.packages('doMC'); library(doMC)} else {library(doMC)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(stringr)){install.packages('stringr'); library(stringr)} else {library(stringr)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(corrplot)){install.packages('corrplot'); library(corrplot)} else {library(corrplot)})
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(factoextra)){install.packages('factoextra'); library(factoextra)} else {library(factoextra)})
suppressMessages(if(!require(leaflet)){install.packages('leaflet'); library(leaflet)} else {library(leaflet)})
suppressMessages(if(!require(Rtsne)){install.packages('Rtsne'); library(Rtsne)} else {library(Rtsne)})
suppressMessages(if(!require(readr)){install.packages('readr'); library(readr)} else {library(readr)})
suppressMessages(if(!require(dbscan)){install.packages('dbscan'); library(dbscan)} else {library(dbscan)})
suppressMessages(if(!require(zoom)){install.packages('zoom'); library(zoom)} else {library(zoom)})


crop <- "Bean"
occ_data <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Bean/database/occ_data_sum.csv")
crop_cycle_id <- function(crop, occ_data){
  
  crop_list <- tolower(list.files(path = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data", full.names = F))
  ggcmi_lst <- c("Pulses", NA, "Barley", NA, rep("Pulses", 3), NA, "Millet", rep("Pulses", 2), "Millet", "Pulses", "Potatoes", "Rice", "Sorghum", "Sunflower", NA, "Wheat")
  if(crop %in% crop_list){
    crop_ggcmi <- ggcmi_lst[which(crop_list == crop)]
  }
  
  # Planting dates
  planting_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/", crop_ggcmi, "_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
  planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
  
  # Harversting dates
  harvest_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/", crop_ggcmi, "_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
  harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
  
  planting <- raster::extract(x = planting_rf_ggcmi, y = occ_data[,c("lon", "lat")])
  harvest <- raster::extract(x = harvest_rf_ggcmi, y = occ_data[,c("lon", "lat")])
  occ_data <- cbind(occ_data, planting, harvest); rm(planting, harvest, planting_rf_ggcmi, harvest_rf_ggcmi)
  occ_data$planting[which(occ_data$planting == "-99")] <- NA
  occ_data$harvest[which(occ_data$harvest == "-99")] <- NA
   
  if(!file.exists(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/precense_data/",crop_ggcmi,"/plots/summary_of_missing_data_/",crop_ggcmi," .png"))){
    ## Proporcion de NA 
   total <-  length(occ_data$planting)
   x <- sum(is.na(occ_data$planting) & !is.na(occ_data$harvest))/total
   y <- sum(!is.na(occ_data$planting) & is.na(occ_data$harvest))/total
   z <- sum(is.na(occ_data$planting) & is.na(occ_data$harvest))/total

  df <- data.frame(trt = c('Planting','Harvest','Planting & Harvest'), outcome= c(x,y,z))
  gg <- ggplot(df, aes(trt, outcome)) +geom_col()
  gg <- gg + labs(title="Summary of Missing Data", y="Proportion of Na", x= "Variables")
  gg
  rm(x,y,z,total,df,gg)
  } 
 
  # Assuming complete data we will continue with the process
  occ_data <- occ_data[complete.cases(occ_data),]; rownames(occ_data) <- 1:nrow(occ_data)
  # Determine length of cycle
  
  occ_data$cycle_length <- abs(occ_data$harvest - occ_data$planting)
  occ_data$condition  <- NA
  occ_data$condition[which((occ_data$harvest - occ_data$planting)> 0)] <- "One year"
  occ_data$condition[which((occ_data$harvest - occ_data$planting)< 0)] <- "Two years"
  

  occ_data1 <- occ_data[,c("lon","lat","planting","harvest","cycle_length")]
  write.csv(occ_data, file("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/database/occ_data_full.csv"))
}
