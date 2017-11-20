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


###########TEST 
####PREC 
america <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_america.rds")
africa <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_africa.rds")
asia <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_asia.rds")
europa <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_europa.rds")
oceania <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered_oceania.rds")

#2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
##########
#  LINUX #
#########
prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds'))
fil <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Bean/database/occ_data_complete_case.csv")
fil <- data.frame(cellId = fil$cellID, lon = fil$lon , lat = fil$lat)

cellID <- c (cellFromXY(base, cbind(x=fil$lon, y= fil$lat)))
y <- match(cellID, prec$cellID)
y <- na.omit(y)
prec1 <- prec[y,]


prec1 <- filter(prec,  prec$cellID  %in% fil$cellID)

# Planting dates
planting_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(planting_rf_ggcmi)
prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(harvest_rf_ggcmi)
prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")

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
}, mc.cores = 20, mc.preschedule = F))


tabla <- do.call(rbind,prec_indexes)

tabla <- rbind(prec_indexes)
##=-=-=-=-=-=-=-=-=-=-##
## conjunto de indices##   Estandarizados
##=-=-=-=-=-=-=-=-=-=-##
tabla <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind.rds"))
tabla <- filter(tabla, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95"  )
tabla$Year <- as.numeric(tabla$Year)
tabla$Value <- as.numeric(tabla$Value)
tabla$cellID <- as.numeric(tabla$cellID)
vars <- unique(tabla$Variable)

spread_tables <- lapply(1:length(vars), function(i){
  df <- tabla %>% filter(Variable == vars[i])
  df <- df %>% dplyr::group_by(Variable) %>% tidyr::spread(Year, Value)
  df$Variable <- NULL
  names(df)[2:ncol(df)] <- paste0(vars[i], "-", names(df)[2:ncol(df)])
  df2 <- apply(X=df[,-1], MARGIN = 1, scale)
  df2 <- t(df2)
  df2<- data.frame(cellID=df[,1],df2)
  names(df2)[2:ncol(df2)] <-  names(df)[2:ncol(df)]
  return(df2)
})
tabla2 <- Reduce(function(...) merge(..., by = "cellID", all.x= T), spread_tables)

indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada CIAT cellID 
cellIdentified <- cellFromXY(base, xy= c(-76.356751, 3.504854)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)
tabla2[is.na(tabla2)] <- 0

dissMat <- matrix(data = NA, nrow = nrow(tabla2), ncol = 3)
system.time (for(i in 1:nrow(tabla2)){
  dissMat[i, 1] <- TSclust::diss.EUCL(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
  dissMat[i, 2] <- TSclust::diss.DTWARP(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
  dissMat[i, 3] <- TSclust::diss.CORT(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]), k = 2, deltamethod = "Euclid")
}); rm(i)






suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
?suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
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

prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered_america.rds'))
fil <- read.csv(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/occ_data_complete_case.csv"))
fil <- data.frame(cellId = fil$cellID, lon = fil$lon , lat = fil$lat)

cellID <- c (cellFromXY(base, cbind(x=fil$lon, y= fil$lat)))
y <- match(prec$cellID,cellID)


prec2 <- prec[y,]
y <- prec2[complete.cases(prec2$cellID),]



############################################################################
tabla <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/tabla_rbind.rds"))
tabla <- filter(tabla, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95"  )
tabla$Year <- as.numeric(tabla$Year)
tabla$Value <- as.numeric(tabla$Value)
tabla$cellID <- as.numeric(tabla$cellID)
vars <- unique(tabla$Variable)

spread_tables <- lapply(1:length(vars), function(i){
  df <- tabla %>% filter(Variable == vars[i])
  df <- df %>% dplyr::group_by(Variable) %>% tidyr::spread(Year, Value)
  df$Variable <- NULL
  names(df)[2:ncol(df)] <- paste0(vars[i], "-", names(df)[2:ncol(df)])
  df2 <- apply(X=df[,-1], MARGIN = 1, scale)
  df2 <- t(df2)
  df2<- data.frame(cellID=df[,1],df2)
  names(df2)[2:ncol(df2)] <-  names(df)[2:ncol(df)]
  return(df2)
})
tabla2 <- Reduce(function(...) merge(..., by = "cellID", all.x= T), spread_tables)

indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada CIAT cellID 
cellIdentified <- cellFromXY(base, xy= c(-76.356751, 3.504854)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)
tabla2[is.na(tabla2)] <- 0

cellID <- c(cellFromXY(base, cbind(x=fil$lon, y= fil$lat)))
y <- match(cellID,tabla2$cellID)
y1 <-na.omit(y)
tab <- tabla2[y1,]
tab1 <- tab[!duplicated(tab$cellID), ]
tabla2 <- tab1
####################################
dissMat <- matrix(data = NA, nrow = nrow(tabla2), ncol = 3)
system.time (for(i in 1:nrow(tabla2)){
  dissMat[i, 1] <- TSclust::diss.EUCL(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
  dissMat[i, 2] <- TSclust::diss.DTWARP(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
  dissMat[i, 3] <- TSclust::diss.CORT(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]), k = 2, deltamethod = "Euclid")
}); rm(i)
dissMat <- as.data.frame(dissMat)
#dissMat[i, 3] <- TSclust::diss.FRECHET(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
colnames(dissMat) <- c("Eucl", "DTWarp", "CORT")
require(tidyr)
tdist <- data.frame(cellID = tabla2$cellID, dissMat)
saveRDS(tdist, paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
####################################################################################


countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
america <- countries[countries@data$CONTINENT == "South America"|
                       countries@data$CONTINENT =="North America"|
                       countries@data$UNREG1 =="Central America",]

###################################################################

dist <- readRDS(paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
base[] <- NA

dist2 <- dist
dist2$Eucl <- cut(dist2$Eucl, quantile(dist2$Eucl, prob = seq(0, 0.5, by = .1)))
dist2$Eucl <- as.character(dist2$Eucl)
dist2$Eucl[which(dist2$Eucl == "(0,11.9]")] <- "10%"
dist2$Eucl[which(dist2$Eucl == "(11.9,13.3]")] <- "20%"
dist2$Eucl[which(dist2$Eucl == "(13.3,14]")] <- "30%"
dist2$Eucl[which(dist2$Eucl == "(14,14.4]")] <- "40%"
dist2$Eucl[which(dist2$Eucl == "(14.4,14.7]")] <- "50%"
dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
myPalette <- colorRampPalette(rev(brewer.pal(5, "Spectral")))
ggplot(data = dist2, aes(x = x, y = y, fill = Eucl)) +
  #scale_fill_gradientn(colours = myPalette(5)) +
  # scale_fill_distiller(palette = "Spectral", na.value="white") +
  geom_raster() +
  coord_equal() +
  theme_bw()
# Hasta aqui funciona
#############################################################################################################


dist <- cbind(dist, xyFromCell(object = base, cell = dist$cellID))
dist <- dist %>% select(cellID, x, y, Eucl, DTWarp, CORT)
dist <- dist %>% gather(key = Metric, value = Value, -(cellID:y))


##### Test
d <- density(dissMat$Eucl)
plot(d)
abline(v = quantile(dissMat$Eucl, prob = seq(0, 0.5, by = .1)), col = 4)



# add a color map with 5 colors
col = terrain.colors(6)
plot(base, col=col, breaks=quantile(dissMat$Eucl, prob = seq(0, 0.5, by = .1)), main="", xlim = c(-180, 0))

# One plot
breaks <- quantile(dissMat$DTWarp, prob = seq(0, 0.5, by = .1))
myPalette <- colorRampPalette(rev(brewer.pal(6, "Spectral")))
ggplot(data = dist %>% filter(Metric == "DTWarp"), aes(x = x, y = y, fill = Value)) +
  scale_fill_gradientn(colours = myPalette(6),#c("red","green","blue","black","pink"),
                       breaks = breaks, labels = format(breaks)) +
  # scale_fill_distiller(palette = "Spectral", na.value="white") +
  geom_raster() +
  coord_equal() +
  theme_bw()

# All plots
suppressMessages(library(ggplot2))
plot_func <- function(df, name) {
  ggplot(data = df, aes(x = x, y = y, fill = Value)) +
    geom_raster() +
    coord_equal() +
    theme_bw() +
    scale_fill_continuous(name = name)
}

# plot_func(df = soil_data, name = "Soilcp")
nested_tmp <- dist %>% 
  group_by(Metric) %>% 
  nest() %>% 
  mutate(plots = map2(data, Metric, plot_func)) 

gridExtra::grid.arrange(grobs = nested_tmp$plots)





theme_set(theme_bw())
X1<- gplot(base)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral", na.value="white") +geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Distance Eclidea")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()


X1<- gplot(base)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_colour_gradientn(colours = c("red","white","blue"),
                         breaks = quantile(dissMat$Eucl, prob = seq(0, 0.5, by = .1)), labels = format(quantile(dissMat$Eucl, prob = seq(0, 0.5, by = .1)))) +
  geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Distance Eclidea")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()




###################

dist <- readRDS(paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
base1 <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
results1 <- data.frame (cellID= dist$cellID, DTW= dist$DTW) #  dist DTW
a <- results1[!duplicated(results1$cellID), ]
summary(a)
min <- 13.68
max <-  15.60
a$condition <- NA
a$condition[which(a$DTW <= min)] <- 1 
a$condition[which(a$DTW >= max)] <-2 
a$condition[which((a$DTW > min) & (a$DTW < max) )] <- 3
a<- na.omit(a)
as.factor(a$condition)

base1[] <- NA
base1[][a$cellID] <- a$condition
base1 <- raster::crop(x = base, y = extent(america)) # Crop according to the study region
library(ggplot2)

theme_set(theme_bw())
X2 <- gplot(base1)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral", na.value="white") +geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Distance DTW")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()
X2
##################
dist <- readRDS(paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
base1 <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
results1 <- data.frame (cellID= dist$cellID, CORT= dist$CORT) #  dist DTW
a <- results1[!duplicated(results1$cellID), ]
summary(a)
min <- 13.68
max <-  15.60
a$condition <- NA
a$condition[which(a$CORT <= min)] <- 1 
a$condition[which(a$CORT >= max)] <-2 
a$condition[which((a$CORT > min) & (a$CORT < max) )] <- 3
a<- na.omit(a)
as.factor(a$condition)

base1[] <- NA
base1[][a$cellID] <- a$condition
base1 <- raster::crop(x = base, y = extent(america)) # Crop according to the study region
library(ggplot2)

theme_set(theme_bw())
X3 <- gplot(base1)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral", na.value="white") +geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Index CORT")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()
X3

library(ggplot2)
library(gridExtra)
grid.arrange(X1,X2, X3, ncol=2, nrow=2,top="Similar Locations to CIAT")
grid.arrange(X2, X3, ncol=2, nrow=1,top="Similar Locations to CIAT")


#########
root <- "//dapadfs/Workspace_cluster_9"

mapspam <- raster::stack(paste(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/bean_mapspam.nc", sep = "")) 
monfreda <- raster::brick(paste(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/bean_monfreda.nc", sep = ""), lvar=4)
occ_data <- read.csv(paste(root, "/CWR_pre-breeding/Input_data/presence_data/", str_to_title(crop), "/database/", crop, "_genesys.csv", sep = ""))
crop_info <- list(mapspam, monfreda, occ_data); names(crop_info) <- c("mapspam", "monfreda", "occ_data")
rm(mapspam, monfreda, occ_data)

a <- cbind(a, xyFromCell(object = base, cell = a$cellID))
a <- cbind(a, raster::extract(x = monfreda[[1]], a[,c("x", "y")]))
names(a)[ncol(a)] <- "Harvested_area"

data4hist <- a %>% group_by(condition) %>% summarise(sum(Harvested_area*1000, na.rm=T))
data4hist$condition[1] <- "High similarity"
data4hist$condition[2] <- "Medium similarity"
data4hist$condition[3] <- "Low similarity"
data4hist$condition <- factor(x = data4hist$condition, levels = c("Low similarity", "Medium similarity", "High similarity"))
names(data4hist)[2] <- "Area"
ggplot(data4hist, aes(x = condition, y = Area, fill = condition)) + geom_col() +
  xlab("Similarity") + ylab("Area (ha)") + scale_fill_brewer(palette = "Set1")

# Use a threshold to define crop areas
crop_info$mapspam[which(crop_info$mapspam[] <= 0)] <- NA
crop_info$mapspam[which(crop_info$mapspam[] > 0)] <- 1

crop_info$monfreda[which(crop_info$monfreda[] <= 0)] <- NA
crop_info$monfreda[which(crop_info$monfreda[] > 0)] <- 1


#########################


dist <- readRDS(paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
base1 <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
results1 <- data.frame (cellID= dist$cellID, DTW= dist$DTW) #  dist DTW
a <- results1[!duplicated(results1$cellID), ]
quantile(a$DTW, prob = seq(0, 1, length = 11), type = 5)
min <- 11.94410
max <-  15.08047
a$condition <- NA
a$condition[which(a$DTW <= min)] <- 1 
a$condition[which(a$DTW >= max)] <-2 
a$condition[which((a$DTW > min) & (a$DTW < max) )] <- 3
a<- na.omit(a)
as.factor(a$condition)

base1[] <- NA
base1[][a$cellID] <- a$condition
base1 <- raster::crop(x = base, y = extent(america)) # Crop according to the study region
library(ggplot2)

theme_set(theme_bw())
X2 <- gplot(base1)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral", na.value="white") +geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Distance DTW (Decil 10%)")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()
X2

######################
dist <- readRDS(paste0( root, "/CWR_pre-breeding/Results/input_tables/chirps/tdist_america.rds"))
base1 <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
results1 <- data.frame (cellID= dist$cellID, DTW= dist$DTW) #  dist DTW
a <- results1[!duplicated(results1$cellID), ]
quantile(a$DTW, prob = seq(0, 1, length = 11), type = 5)
quantile(a$DTW, prob = seq(0, 1, by = 0.1))
min <- 11.94445# 13.31153
max <-  15.08046
a$condition <- NA
a$condition[which(a$DTW <= min)] <- 1 
a$condition[which(a$DTW >= max)] <-2 
a$condition[which((a$DTW > min) & (a$DTW < max) )] <- 3
a<- na.omit(a)
as.factor(a$condition)

base1[] <- NA
base1[][a$cellID] <- a$condition
base1 <- raster::crop(x = base, y = extent(america)) # Crop according to the study region
library(ggplot2)

theme_set(theme_bw())
X3<- gplot(base1)+geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral", na.value="white") +geom_polygon(data=america, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  ggtitle("Distance DTW (Decil 20%)")+
  scale_x_continuous(limits = c(-180, 0))+
  coord_equal()
X3


