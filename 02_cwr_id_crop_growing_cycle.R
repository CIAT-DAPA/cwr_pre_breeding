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


crop <- "potato"
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

  write.csv(occ_data,file = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/",crop_ggcmi,"/occ_data.csv", row.names = F)
  
  
  set.seed (1235)
  smplID <- sample(x = 1:nrow(occ_data), size = 2000, replace = F)
  occ_data_smpl <- occ_data[smplID,]; rownames(occ_data_smpl) <- 1:nrow(occ_data_smpl)
  M<- cor(occ_data_smpl[,c("lon","lat", "planting","harvest","cycle_length")],method = "spearman")
  corrplot(M, method = "ellipse")
  res_pca  <- FactoMineR::PCA(occ_data_smpl[,c("lon","lat", "planting","harvest","cycle_length")], graph = F)
  set.seed(1235)
  res_hcpc <- FactoMineR::HCPC(res_pca, nb.clust = -1, graph = F) # Define number of SCENARIOS automatically
  test <- data.frame(cellID = occ_data_smpl[,1], res_hcpc$data.clust, condition = occ_data_smpl[,7])
 
  if(!file.exists(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/precense_data/",crop_ggcmi,"/plots/boxplots_for_set_of_variables_/",crop_ggcmi," .png"))){
  df_cluster <-res_hcpc$data.clust 
  df_cluster$Combination <- rownames(df_cluster)
  df_cluster <- df_cluster %>% gather(Variables, Value, lon:cycle_length)
  gg <- ggplot(df_cluster, aes(x=clust, y=Value)) + geom_boxplot() 
  gg <- gg + facet_wrap(~ Variables, scales='free_y')
  gg <- gg + ylab('Values') + xlab('Cluster')
  gg
  rm(df_cluster, gg)
  }
  
  raster_occ_data<- raster::stack(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/potato/database/potato_monfreda.nc", sep = ""))
  raster_occ_data <- raster_occ_data[[1]]
  raster_occ_data[which(raster_occ_data[] > 0)] <- NA
  raster_occ_data[which(raster_occ_data[] <= 0)] <- NA
  #points(x=res_hcpc$data.clust$lon,y= res_hcpc$data.clust$lat,col = res_hcpc$data.clust$clust, pch = 20)
  raster_occ_data[test$cellID]<-as.numeric(test$clust)
  plot(raster_occ_data)
 
  
  r <- raster_occ_data
  pal <- colorNumeric(c("red", "blue", "yellow", "pink"), values(r),
                      na.color = "transparent")
  leaflet() %>% addTiles() %>%
    addRasterImage(r, colors = pal, opacity = 0.8) %>%
    addLegend(pal = pal, values = values(r),
              title = "Surface temp")
  
   # Prueba numero de muestra
   lista <- c()
    n <- c(100,200,300,500,800,1000,1500,2000,2500,3000,3500,4000,5000,6000,7000)
   for (i in 1:length(n)){
   set.seed (1235)
   smplID <- sample(x = 1:nrow(occ_data), size = n[i] , replace = F)
   occ_data_smpl <- occ_data[smplID,]; rownames(occ_data_smpl) <- 1:nrow(occ_data_smpl)
   res_pca  <- FactoMineR::PCA(occ_data_smpl[,c("lon","lat", "planting","harvest","cycle_length")], graph = F)
   set.seed(1235)
    res_hcpc <- FactoMineR::HCPC(res_pca, nb.clust = -1, graph = F) # Define number of SCENARIOS automatically
   test <- data.frame(cellID = occ_data_smpl[,1], res_hcpc$data.clust, condition = occ_data_smpl[,7])
   lista[i] <- length(unique(test$clust))
   }
   
  if(!file.exists(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/precense_data/",crop_ggcmi,"/plots/determine_sample_size_cluster_/",crop_ggcmi," .png"))){
  gg <-  qplot( x=n ,y=lista,xlab= "Sample Size", ylab= "Number of clusters")
  gg <- gg + ggtitle("Determine sample size to do PCA & Cluster analysis")
  gg <- gg + theme(axis.title=element_text(face="bold.italic",size="12", color="brown"), legend.position="top")
  gg
  rm(gg)
  }

set.seed(500)
semi <-c(runif(500,min=1,max=1000000))
semi <- round(semi)   

tablas <- lapply(1:length(semi), function(i){ # 1:length(lista_tmin)
  cat(paste("Procesando semilla:", i, "\n", sep = ""))
  set.seed (semi[i])
  smplID <- sample(x = 1:nrow(occ_data), size = 2000, replace = F)
  occ_data_smpl <- occ_data[smplID,]; rownames(occ_data_smpl) <- 1:nrow(occ_data_smpl)
  M<- cor(occ_data_smpl[,c("lon","lat", "planting","harvest","cycle_length")],method = "spearman")
  res_pca  <- FactoMineR::PCA(occ_data_smpl[,c("lon","lat", "planting","harvest","cycle_length")], graph = F)
  set.seed(semi[i])
  res_hcpc <- FactoMineR::HCPC(res_pca, nb.clust = -1, graph = F) # Define number of SCENARIOS automatically
  test <- data.frame(cellID = occ_data_smpl[,1], res_hcpc$data.clust, condition = occ_data_smpl[,7])
   return(test)
  removeTmpFiles(h=0)
})

#points(x=res_hcpc$data.clust$lon,y= res_hcpc$data.clust$lat,col = res_hcpc$data.clust$clust, pch = 20)
raster_occ_data[test$cellID]<-as.numeric(test$clust)
plot(raster_occ_data)

rast <- lapply(1:40, function(i){
cat(paste("Procesando semilla:", i, "\n", sep = ""))
raster_occ_data[tablas[[i]]$cellID]<-as.numeric(tablas[[i]]$clust)
return(raster_occ_data)
removeTmpFiles(h=0)
})  

rastsemi <- raster::brick(rast) 


r <- rastsemi
pal <- colorNumeric(c("red", "blue", "yellow", "pink"), values(r),
                    na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(r, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r),
            title = "Surface temp")



}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 

###cc_scenario  line 475-510
#m <- leaflet() %>%
  #addTiles() %>%  # Add default OpenStreetMap map tiles
  #addMarkers(lng= res_hcpc$data.clust$lon, lat=res_hcpc$data.clust$lat, popup="The birthplace of R")
#m  # Print the map



### Sotelo, shape, validation, 

