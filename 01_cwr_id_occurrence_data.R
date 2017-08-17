# CWR pre-breeding characterising testing environments: Identify crop areas
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
##   Extraction of raster information   ##
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

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
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(rasterVis)){install.packages('rasterVis'); library(rasterVis)} else {library(rasterVis)})
suppressMessages(if(!require(stringr)){install.packages('stringr'); library(stringr)} else {library(stringr)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(factoextra)){install.packages('factoextra'); library(factoextra)} else {library(factoextra)})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#   Load Monfread and Mapspam   #
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

crop_area_id <- function(crop){
  
  OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
  if(OSys == "Linux"){ root <- "/mnt/workspace_cluster_9" } else {
    if(OSys == "Windows"){ root <- "//dapadfs/Workspace_cluster_9" }
  }; rm(OSys)
  
  # Load area information by crop
  mapspam <- raster::stack(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/database/", crop, "_mapspam.nc", sep = ""))
  monfreda <- raster::stack(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/database/",crop, "_monfreda.nc", sep = ""))
  occ_data <- read.csv(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/database/", crop, "_genesys.csv", sep = ""))
  crop_info <- list(mapspam, monfreda, occ_data); names(crop_info) <- c("mapspam", "monfreda", "occ_data")
  rm(mapspam, monfreda, occ_data)
  
  # Use a threshold to define crop areas
  crop_info$mapspam[which(crop_info$mapspam[] <= 0)] <- NA
  crop_info$mapspam[which(crop_info$mapspam[] > 0)] <- 1
  
  crop_info$monfreda[which(crop_info$monfreda[] <= 0)] <- NA
  crop_info$monfreda[which(crop_info$monfreda[] > 0)] <- 1
  
  # Plotting maps of occurrence data
  tmpStack <- raster::stack(crop_info$mapspam,
                            crop_info$monfreda)
  rasterSum <- sum(tmpStack, na.rm = T)
  rasterSum[rasterSum[] == 0] <- NA
  tmpStack[[3]] <- rasterSum; rm(rasterSum)
  names(tmpStack) <- c("MapSPAM", "Monfreda", "MapSPAM_Monfreda")
  
  # Make the plot
  if(!file.exists(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots/occurrence_data_", crop, "_maps.png", sep = ""))){
    detach(package:factoextra)
    detach(package:ggplot2)
    
    # Load shapefile worldwide
    data(wrld_simpl)
    
    trellis.device(device="png", filename= paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots/occurrence_data_", crop, "_maps.png", sep = ""), width = 3000, height = 2000, units = "px", res = 300)
    mapTheme <- rasterTheme(region = brewer.pal(3,"BrBG")) ## gusta 2
    p <- levelplot(tmpStack, margin = F, par.settings = mapTheme) + layer(sp.polygons(wrld_simpl, lwd = 0.1, col = 'gray')) #c('palegreen', 'indianred1')
    print(p)
    dev.off()
    rm(p, wrld_simpl)
    
  }
  
  # Determine agreement of each raster file with Genesys database
  if(!file.exists(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots/occurrence_data_", crop, "_agreement.png", sep = ""))){
    agreeMapspam <- raster::extract(x = tmpStack[[1]], y = crop_info$occ_data[,c("longitude","latitude")])
    agreeMapspam <- sum(agreeMapspam, na.rm = T)/nrow(crop_info$occ_data)
    agreeMonfreda <- raster::extract(x = tmpStack[[2]],y = crop_info$occ_data[,c("longitude","latitude")])
    agreeMonfreda <- sum(agreeMonfreda, na.rm = T)/nrow(crop_info$occ_data)
    agreeSum <- raster::extract(x = tmpStack[[3]],y = crop_info$occ_data[,c("longitude","latitude")])
    agreeSum <- length(agreeSum[which(agreeSum == 2)])/nrow(crop_info$occ_data)
    
    agree <- data.frame(Dataset = factor(x = c("MapSPAM", "Monfreda", "MapSPAM & Monfreda"), levels = c("MapSPAM", "Monfreda", "MapSPAM & Monfreda")), Percentage = c(agreeMapspam, agreeMonfreda, agreeSum))
    rm(agreeMapspam, agreeMonfreda, agreeSum)
    
    if(!dir.exists(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots", sep = ""))){
      dir.create(path = paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots", sep = ""), recursive = T)
    }
    gg <- agree %>% ggplot(aes(x=Dataset, y=Percentage * 100)) + geom_bar(stat="identity", fill="steelblue") + theme_bw()
    gg <- gg + ylab("Percentage (%)") + ggtitle(label = "Agreement between Genesys and raster sources")
    gg <- gg + scale_y_continuous(limits = c(0, 100))
    ggsave(filename = paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/plots/occurrence_data_", crop, "_agreement.png", sep = ""), plot = gg, width = 5, height = 5.5, units = "in")
    rm(gg, agree, tmpStack)
  }
  
  # Just for now, we will use Monfreda (This is a test) Modeling tools are required
  tmpStack$MapSPAM_Monfreda[which(tmpStack$MapSPAM_Monfreda[]==2)] <- 1
  data_matrix <- as.data.frame(rasterToPoints(tmpStack$MapSPAM_Monfreda))
  cellID <- cellFromXY(object =tmpStack$MapSPAM_Monfreda, xy = data_matrix[,c("x", "y")])
  data_matrix <- cbind(cellID, data_matrix[,c("x", "y")]); colnames(data_matrix) <- c("cellID", "lon", "lat")
  occ_data <- data_matrix
  write.csv(occ_data, paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/database/occ_data_sum.csv", sep = ""), row.names = F)
  
  return(data_matrix)
  
}
