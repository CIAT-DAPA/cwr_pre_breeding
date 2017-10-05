# CWR pre-breeding characterising testing environments: Index calculation for current data
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

##=-=-=-=-=-=-=-=-=-=-##
## Correr en Windows  ##
##=-=-=-=-=-=-=-=-=-=-##

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

# Verify coordinates within Central America, Caribbean and Colombia
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
central_america_colombia <- countries[countries@data$COUNTRY == "Colombia" |
                                        countries@data$COUNTRY == "Venezuela" |
                                        countries@data$COUNTRY == "Ecuador"|
                                        countries@data$COUNTRY == "Peru" |
                                        countries@data$COUNTRY == "Bolivia" |
                                        countries@data$COUNTRY == "Argentina" |
                                        countries@data$COUNTRY == "Chile" |
                                        countries@data$COUNTRY == "Brazil" |
                                        countries@data$UNREG1 == "Central America" |
                                        countries@data$UNREG1 == "Caribbean",]
# Filter coordinates within Central America, Caribbean and Colombia for all variables
varList <- c("prec", "srad", "tmax", "tmin")
for(i in 1:length(varList)){
  eval(parse(text = paste0("over_res <- sp::over(SpatialPoints(coords = data.frame(lon = ", varList[i], "$lon, lat = ", varList[i], "$lat), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')), as(central_america_colombia, 'SpatialPolygons'))")))
  eval(parse(text = paste0(varList[i], "$bean_coordinates <- over_res; rm(over_res)")))
  eval(parse(text = paste0(varList[i], " <- ", varList[i], "[!is.na(", varList[i], "$bean_coordinates),]; rownames(", varList[i], ") <- 1:nrow(", varList[i], ")")))
}; rm(i, varList)

saveRDS(prec, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered.rds"))
saveRDS(srad, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_filtered.rds"))
saveRDS(tmax, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_filtered.rds"))
saveRDS(tmin, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_filtered.rds"))

###########################
#Cargar bases simultaneas#
##########################

##=-=-=-=-=-=-=-=-=-=-##
##  Correr en LINUX   ##
##=-=-=-=-=-=-=-=-=-=-##

system.time(prec <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_filtered.rds")))
prec$bean_coordinates <- NULL

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

##=-=-=-=-=-=-=-=-=-=-##
## Correr en Windows  ##  POINT
##=-=-=-=-=-=-=-=-=-=-##
library(Rcpp)
library(RcppArmadillo)
library(raster)
library(dplyr)
library(rgdal)

##=-=-=-=-=-=-=-=-=-=-##
# Indices uno por uno #
##=-=-=-=-=-=-=-=-=-=-##

# OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
# if(OSys == "Linux"){
#   root <- "/mnt/workspace_cluster_9"
#   base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
# } else {
#   if(OSys == "Windows"){
#     root <- "//dapadfs/Workspace_cluster_9"
#     base <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/AgMerra_template.RDS"))
#   }
# }; rm(OSys)

tabla <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind.rds"))
require(dplyr)
a <- filter(tabla, Variable  == "Error in summarise_impl(.data, dots) : \n  Evaluation error: only 0's may be mixed with negative subscripts.\n")
b <- filter(tabla, Variable  == "Error in xj[i] : only 0's may be mixed with negative subscripts\n")
rm(a,b)
index <- c("TOTRAIN","CDD","P5D","P_95")
sourceCpp("fastPdist.cpp") # Pairwise Euclidean distance function in C++
tabla <- filter(tabla, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95"  )
tabla$Year <- as.numeric(tabla$Year)
tabla$Value <- as.numeric(tabla$Value)

countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
central_america_colombia <- countries[countries@data$COUNTRY == "Colombia" |
                                        countries@data$COUNTRY == "Venezuela" |
                                        countries@data$COUNTRY == "Ecuador"|
                                        countries@data$COUNTRY == "Peru" |
                                        countries@data$COUNTRY == "Bolivia" |
                                        countries@data$COUNTRY == "Argentina" |
                                        countries@data$COUNTRY == "Chile" |
                                        countries@data$COUNTRY == "Brazil" |
                                        countries@data$UNREG1 == "Central America"|
                                        countries@data$UNREG1 == "Caribbean",]


ciclo <- lapply(1:length(index),function(i){
  cat(paste0("procesando  index", i,"\n")) 
  indexes <- tabla%>% filter( Variable == index[i])
  indexes$Variable <- NULL
  require(tidyr)
  indexes <- indexes %>% spread(key = Year, value = Value)
  rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
  ###coordenada CIAT cellID 
  cellIdentified <- cellFromXY(base, xy= c(-76.356751, 3.504854)) ## Coordenada a buscar
  cor<- which(rownames(indexes) == cellIdentified)
  pdist_result <- fastPdist2(Ar = as.matrix(indexes), Br = as.matrix(indexes[cor,])) # Here we select the pre-breeding site
  results <- data.frame(cellID = rownames(indexes), pdist = as.numeric(pdist_result))
  base[] <- NA
  base[][results$cellID %>% as.character %>% as.numeric] <- results$pdist
  
  base <- raster::crop(x = base, y = extent(central_america_colombia)) # Crop according to the study region
  return(base)
})
##NOTA i= 1,2,3,4
#1. TOTRAIN
#2. CDD
#3. P5D
#4. P_95

require(levelplot)
levelplot(ciclo[[i]], par.settings = RdBuTheme)
### grafico ggplot
require(ggplot2)
theme_set(theme_bw())
gplot(ciclo[[2]])+ geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral") +geom_polygon(data=central_america_colombia, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  coord_equal()

##=-=-=-=-=-=-=-=-=-=-##
## conjunto de indices##   Sin estandarizar 
##=-=-=-=-=-=-=-=-=-=-##

tabla <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind.rds"))
tabla <- filter(tabla, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95"  )
tabla$Year <- as.numeric(tabla$Year)
tabla$Value <- as.numeric(tabla$Value)
vars <- unique(tabla$Variable)
spread_tables <- lapply(1:length(vars), function(i){
  df <- tabla %>% filter(Variable == vars[i])
  df <- df %>% dplyr::group_by(Variable) %>% tidyr::spread(Year, Value)
  df$Variable <- NULL
  names(df)[2:ncol(df)] <- paste0(vars[i], "-", names(df)[2:ncol(df)])
  return(df)
})
tabla2 <- Reduce(function(...) merge(..., by = "cellID", all.x = T), spread_tables)
indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada CIAT cellID 
cellIdentified <- cellFromXY(base, xy= c(-76.356751, 3.504854)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)
sourceCpp("fastPdist.cpp") # Pairwise Euclidean distance function in C++
pdist_result <- fastPdist2(Ar = as.matrix(indexes), Br = as.matrix(indexes[cor,])) # Here we select the pre-breeding site
results <- data.frame(cellID = rownames(indexes), pdist = as.numeric(pdist_result))
base[] <- NA
base[][results$cellID %>% as.character %>% as.numeric] <- results$pdist
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
central_america_colombia <- countries[countries@data$COUNTRY == "Colombia"|
                                        countries@data$COUNTRY == "Venezuela"|
                                        countries@data$COUNTRY == "Ecuador"|
                                        countries@data$COUNTRY == "Peru"|
                                        countries@data$COUNTRY == "Bolivia"|
                                        countries@data$COUNTRY == "Argentina"|
                                        countries@data$COUNTRY == "Chile"|
                                        countries@data$COUNTRY == "Brazil"|
                                        countries@data$UNREG1 == "Central America"|
                                        countries@data$UNREG1 == "Caribbean",]
base <- raster::crop(x = base, y = extent(central_america_colombia)) # Crop according to the study region
plot(base)


library(ggplot2)
theme_set(theme_bw())
gplot(base)+ geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral") +geom_polygon(data=central_america_colombia, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  coord_equal()
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
sourceCpp("fastPdist.cpp") # Pairwise Euclidean distance function in C++
pdist_result <- fastPdist2(Ar = as.matrix(indexes), Br = as.matrix(indexes[cor,])) # Here we select the pre-breeding site
results <- data.frame(cellID = rownames(indexes), pdist = as.numeric(pdist_result))
base[] <- NA
base[][results$cellID %>% as.character %>% as.numeric] <- results$pdist
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
central_america_colombia <- countries[countries@data$COUNTRY == "Colombia"|
                                        countries@data$COUNTRY == "Venezuela"|
                                        countries@data$COUNTRY == "Ecuador"|
                                        countries@data$COUNTRY == "Peru"|
                                        countries@data$COUNTRY == "Bolivia"|
                                        countries@data$COUNTRY == "Argentina"|
                                        countries@data$COUNTRY == "Chile"|
                                        countries@data$COUNTRY == "Brazil"|
                                        countries@data$UNREG1 == "Central America"|
                                        countries@data$UNREG1 == "Caribbean",]
base <- raster::crop(x = base, y = extent(central_america_colombia)) # Crop according to the study region
plot(base)


library(ggplot2)
theme_set(theme_bw())
gplot(base)+ geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_distiller(palette = "Spectral") +geom_polygon(data=central_america_colombia, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.5)+
  coord_equal()

