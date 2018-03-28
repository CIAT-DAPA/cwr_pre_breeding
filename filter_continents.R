# CWR pre-breeding characterising testing environments: Index calculation for current data
# Authors: H. Achicanoy & B. Mora
# CIAT, 2017

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
#1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

# Precipitation data (from CHIRPS)
prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/table_final.rds'))


# Solar radiation, tmin and tmax data (from AgMerra)
agList <- c("srad", "tmax", "tmin")
for(i in 1:length(agList)){
  
  cat(paste0("\n\n==================== Processing: ", agList[i], "\n"))
  
  eval(parse(text = paste(agList[i], ' <- readRDS(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_", agList[i], "/table_final.rds", sep = ""))', sep = '')))
  eval(parse(text = paste('cellID <- cellFromXY(base, xy = ', agList[i], '[,1:2])', sep = '')))
  eval(parse(text = paste(agList[i], ' <- data.frame(cellID = cellID, ', agList[i], ')', sep = '')))
  eval(parse(text = paste('names(', agList[i], ')[2:3] <- c("lon", "lat")', sep = '')))
  eval(parse(text = paste('names(', agList[i], ')[4:ncol(', agList[i], ')] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(', agList[i], ')[4:ncol(', agList[i], ')]), fixed = T))', sep = '')))
  eval(parse(text = paste(agList[i], ' <- ', agList[i], '[which(', agList[i], '$cellID %in% base::intersect(prec$cellID, ', agList[i], '$cellID)),]; rownames(', agList[i], ') <- 1:nrow(', agList[i], ')', sep = '')))
  rm(cellID)
}; rm(i)

rm(base, agList)
cellID <- srad$cellID
y <- match(cellID, srad$cellID)
y <- na.omit(y)
prec<- prec[y,]



#####################
#    cONTINENTES    #
#####################

# Verify coordinates within Continents 
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

america <- countries[countries@data$CONTINENT == "South America"|
                       countries@data$CONTINENT =="North America"|
                       countries@data$UNREG1 =="Central America",]

africa <- countries[countries@data$CONTINENT == "Africa",]

asia <- countries[countries@data$CONTINENT == "Asia",]
oceania <- countries[countries@data$CONTINENT == "Oceania",]
europa <- countries[countries@data$CONTINENT == "Europe",]

######################
#      AMERICA      #    
######################  

precAm <- prec
sradAm <- srad
tmaxAm <- tmax
tminAm <- tmin 
varList <- c("precAm", "sradAm", "tmaxAm", "tminAm")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(america, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)

saveRDS(precAm, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds"))
saveRDS(sradAm, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_america.rds"))
saveRDS(tmaxAm, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_america.rds"))
saveRDS(tminAm, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_america.rds"))
rm (precAm,sradAm, tmaxAm,tminAm)


######################
#      AFRICA        #    
###################### 
precAf <- prec
sradAf <- srad
tmaxAf<- tmax
tminAf<- tmin    
varList <- c("precAf", "sradAf", "tmaxAf", "tminAf")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(africa, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)

saveRDS(precAf, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_africa.rds"))
saveRDS(sradAf, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_africa.rds"))
saveRDS(tmaxAf, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_africa.rds"))
saveRDS(tminAf, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_africa.rds"))
rm (precAf,sradAf,tmaxAf,tminAf)
######################
#      EUROPA        #    
###################### 
precE <- prec
sradE <- srad
tmaxE <- tmax
tminE <- tmin  
varList <- c("precE", "sradE", "tmaxE", "tminE")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(europa, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)
saveRDS(precE, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_europa.rds"))
saveRDS(sradE, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_europa.rds"))
saveRDS(tmaxE, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_europa.rds"))
saveRDS(tminE, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_europa.rds"))
rm(precE,sradE,tmaxE,tminE)
######################
#        ASIA        #    
######################

precAs <- prec
sradAs <- srad
tmaxAs<- tmax
tminAs <- tmin    
varList <- c("precAs", "sradAs", "tmaxAs", "tminAs")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(asia, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)
saveRDS(precAs, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_asia.rds"))
saveRDS(sradAs, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_asia.rds"))
saveRDS(tmaxAs, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_asia.rds"))
saveRDS(tminAs, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_asia.rds"))
rm(precAs,sradAs,tmaxAs,tminAs)

######################
#     Oceania        #    
######################

precO <- prec
sradO <- srad
tmaxO <- tmax
tminO <- tmin    
varList <- c("precO", "sradO", "tmaxO", "tminO")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(oceania, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)
saveRDS(precO, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_oceania.rds"), compress = TRUE)
saveRDS(sradO, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_oceania.rds"))
saveRDS(tmaxO, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered_oceania.rds"))
saveRDS(tminO, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered_oceania.rds"))
rm(precO,sradO,tmaxO,tminO)


