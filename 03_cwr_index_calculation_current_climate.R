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
suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(FactoClass)){install.packages('FactoClass'); library(FactoClass)} else {library(FactoClass)})
suppressMessages(if(!require(ade4)){install.packages('ade4'); library(ade4)} else {library(ade4)})
suppressMessages(if(!require(xtable)){install.packages('xtable'); library(xtable)} else {library(xtable)})
suppressMessages(if(!require(ggdendro)){install.packages('ggdendro'); library(ggdendro)} else {library(ggdendro)})
suppressMessages(if(!require(compiler)){install.packages('compiler'); library(compiler)} else {library(compiler)})
suppressMessages(if(!require(ggthemes)){install.packages('ggthemes'); library(ggthemes)} else {library(ggthemes)})

OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
if(OSys == "Linux"){
  root <- "/mnt/workspace_cluster_9"
  base <- raster::stack("/mnt/data_cluster_5/cropdata/agmerra/daily/nc-files/srad_daily_ts_agmerra_1980_2010.nc")
} else {
  if(OSys == "Windows"){
    root <- "//dapadfs/Workspace_cluster_9"
    base <- raster::stack("//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/srad_daily_ts_agmerra_1980_2010.nc")
  }
}; rm(OSys)

base <- base[[1]]
base <- rotate(base)

# Precipitation data (from CHIRPS)
prec <- readRDS(paste0(root, '/CWR_pre-breeding/Results/input_tables/chirps/table_final.rds'))
cellID <- cellFromXY(base, xy = prec[,1:2])
prec <- data.frame(cellID = cellID, prec); rm(cellID)
names(prec)[2:3] <- c("lon", "lat")
names(prec)[4:ncol(prec)] <- as.character(gsub(pattern = "chirps.v2.0.", replacement = "", x = names(prec)[4:ncol(prec)]))
names(prec)[4:ncol(prec)] <- gsub(pattern = "\\.", replacement = "-", x = names(prec)[4:ncol(prec)])

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
central_america_colombia <- countries[countries@data$COUNTRY == "Colombia" | countries@data$UNREG1 == "Central America" | countries@data$UNREG1 == "Caribbean",]
over_res <- sp::over(SpatialPoints(coords = data.frame(lon = prec$lon, lat = prec$lat), proj4string = CRS(projargs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")), as(central_america_colombia, "SpatialPolygons"))
prec$bean_coordinates <- over_res; rm(over_res)
prec <- prec[!is.na(prec$bean_coordinates),]; rownames(prec) <- 1:nrow(prec)

saveRDS(prec, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/prec_central_america_colombia.rds"))
saveRDS(srad, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_new.rds"))
saveRDS(tmax, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_new.rds"))
saveRDS(tmin, paste0(root, "/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_new.rds"))

rm(eli,eli1,prec,ggcmi_data)






prec <- readRDS("D:/prec_CACol.rds")
prec$bean_coordinates <- NULL

# Planting dates
planting_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(planting_rf_ggcmi)
prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")]); rm(harvest_rf_ggcmi)
prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")

prec[1:5,(ncol(prec)-5):ncol(prec)]

table(prec$Duration)

saveRDS(object = prec, file = "D:/precipitation.rds")

prec <- readRDS("precipitation.rds")

TEST <- lapply(1:nrow(prec), function(i){
  
  duration <- prec$Duration[i]
  start <- prec$Planting[i]
  end <- prec$Harvest[i]
  
  time.serie <- prec[i, 1:(ncol(prec)-3)]
  
  if(duration == "One year"){
    
    calculations <- function(time.serie, start, end){
      
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
      
      return(data.frame(cellID = unique(X$cellID), rbind(totrain, cdd, p5d, p_95)))
      
    }
    results <- calculations(time.serie = time.serie, start = start, end = end)
    
  }
  
  return(results)
  
})

TEST <- do.call(rbind, TEST)

TEST %>% ggplot() + theme_tufte() + geom_line(aes(x = Year, y = Value, group = cellID, colour = Variable)) + facet_wrap(~Variable)



















prec <- prec[which(prec$cellID %in% srad$cellID),]; rownames(prec) <- 1:nrow(prec)

prec1 <- intersect(prec[,c("lon","lat")] , df[,c("x","y")])
mapspam   <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/database/potato_mapspam.nc")
monfreda  <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/presence_data/Potato/database/potato_monfreda.nc")

#mapspam
mapspam[which(mapspam[] <= 0)] <- NA
mapspam[which(mapspam[] > 0)] <- 1

#monfreda
monfreda[which(monfreda[] <= 0)] <- NA
monfreda[which(monfreda[] > 0)] <- 1

##superpongo los rasters
tmpStack <- raster::stack(mapspam,monfreda)
# sumo los rasters
rasterSum <- sum(tmpStack, na.rm = T)
rasterSum[rasterSum[] == 0] <- NA
rasterSum[rasterSum[] == 2] <- 1


prec2 <- prec
prec2$cropID <- raster::extract(x = rasterSum, y = prec2[,c("lon", "lat")])
prec2.1 <- prec2[which(!is.na(prec2$cropID)),]
plot(rasterSum)


# Planting dates
planting_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harversting dates
harvest_rf_ggcmi <- raster::brick(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Potatoes_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]


prec2.1$planting_day <- raster::extract(x = planting_rf_ggcmi, y = prec2.1[,c("lon", "lat")])
tail(prec2.1[(10955:10961)])
prec2.1$harvest_day <- raster::extract(x = harvest_rf_ggcmi, y = prec2.1[,c("lon", "lat")])
prec2.1$cycle_length <- abs(prec2.1$harvest_day - prec2.1$planting_day)

rm(coor_prec,coor_prec1, mtcars, first, second,prec2,df,mapspam,monfreda,rasterSum,tmpStack,b )
head(prec2.1[,c(1:3,(ncol(prec2.1)-5):ncol(prec2.1))])

prec2.1$condition  <- NA
prec2.1$condition[which((prec2.1$harvest - prec2.1$planting)> 0)] <- "One year"
prec2.1$condition[which((prec2.1$harvest - prec2.1$planting)< 0)] <- "Two years"

hist(prec2.1$cycle_length)
##### important
prec2.1$planting_day[which(prec2.1$planting_day == "-99")] <- NA
prec2.1$harvest_day[which(prec2.1$harvest_day == "-99")] <- NA
prec2.1 <- prec2.1[which(!is.na(prec2.1$planting_day)),]
prec2.1 <- prec2.1[which(!is.na(prec2.1$harvest_day)),]
##chevere
#prec2.1$planting <- NULL



prec2.1 <- prec2.1[,c("cellID", "lon", "lat", "planting_day", "harvest_day", "cycle_length", "condition", colnames(prec2.1)[8:(ncol(prec2.1)-5)])]
prec2.1_long <- prec2.1 %>% gather(key = Date, value = Value, -(cellID:condition))

two_years <- prec2.1_long[which(prec2.1_long$condition == "Two years"),]
saveRDS(two_years, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/two_years.rds")

one_year <- prec2.1_long[which(prec2.1_long$condition == "One year"),]
# yearList <- 1981:2010
# lapply(yearList, function(x){
#   yearTable <- one_year[which(lubridate::year(one_year$Date) == x),]
#   saveRDS(yearTable, paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/one_year_", x,".rds"))
# })

prec2.1_long <- prec2.1_long %>% arrange(cellID, Date)
prec2.1_long$yday <- lubridate::yday(as.Date(prec2.1_long$Date))
prec2.1_long$Year <- lubridate::year(as.Date(prec2.1_long$Date))
prec2.1_long <- prec2.1_long %>% group_by(cellID, Year) %>% filter(yday >= planting_day & yday <= harvest_day)


c(seq(planting_day:365),seq(1:harvest_day))


cell1_prec2.1_long <- prec2.1_long %>% filter(cellID == 455760)
cell1_prec2.1_long <- cell1_prec2.1_long %>% arrange(cellID, Date)
cell1_prec2.1_long$yday <- lubridate::yday(as.Date(cell1_prec2.1_long$Date))
cell1_prec2.1_long$Year <- lubridate::year(as.Date(cell1_prec2.1_long$Date))
cell1_prec2.1_long <- cell1_prec2.1_long %>% group_by(cellID, Year) %>% filter(yday >= planting_day & yday <= harvest_day)


totrain <- cell1_prec2.1_long %>% dplyr::group_by(cellID, Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
totrain$id_key <- paste(totrain$id_hogar, "-", totrain$Period, sep = "")

gg <- totrain %>% ggplot() + theme_tufte() + geom_line(aes(x = Year, y = TOTRAIN, group = id_key, colour = Period)) + facet_wrap(~municipio) + scale_x_continuous(breaks = 2006:2017) + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "C:/Users/haachicanoy/Documents/CIAT/Asesorias/Lisset_Perez/Agroclimas/new_indices/totalrainFrom2006.png", plot = gg)
rm(tmax,tmin,srad,harvest_rf_ggcmi,planting_rf_ggcmi)



