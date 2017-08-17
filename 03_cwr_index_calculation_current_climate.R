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

OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
if(OSys == "Linux"){ root <- "/mnt/workspace_cluster_9" } else {
  if(OSys == "Windows"){ root <- "//dapadfs/Workspace_cluster_9" }
}; rm(OSys)

base <- raster::stack("//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/srad_daily_ts_agmerra_1980_2010.nc")
base <- base[[1]]
base <- rotate(base)

# Precipitation data (from CHIRPS)

prec <- readRDS('//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/table_final.rds')
cellID <- cellFromXY(base, xy = prec[,1:2])
prec <- data.frame(cellID = cellID, prec); rm(cellID)
names(prec)[2:3] <- c("lon", "lat")
names(prec)[4:ncol(prec)] <- as.character(gsub(pattern = "chirps-v2-0-", replacement = "", x = names(prec)[4:ncol(prec)]))

# Solar radiation, tmin and tmax data (from AgMerra)

agList <- c("srad", "tmax", "tmin")

# lapply(1:length(agList), function(i){
#   
#   cat(paste0("\n\n==================== Processing: ", agList[i], "\n"))
#   
#   eval(parse(text = paste(agList[i], ' <- readRDS(paste("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_", agList[i], "/table_final.rds", sep = ""))', sep = '')))
#   eval(parse(text = paste('cellID <- cellFromXY(base, xy = ', agList[i], '[,1:2])', sep = '')))
#   eval(parse(text = paste(agList[i], ' <- data.frame(cellID = cellID, ', agList[i], ')', sep = '')))
#   eval(parse(text = paste('names(', agList[i], ')[2:3] <- c("lon", "lat")', sep = '')))
#   eval(parse(text = paste('names(', agList[i], ')[4:ncol(', agList[i], ')] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(', agList[i], ')[4:ncol(', agList[i], ')]), fixed = T))', sep = '')))
#   eval(parse(text = paste(agList[i], '<-', agList[i], '[', agList[i], '$CellID %in% intersect(prec$cellID, ', agList[i], '$cellID)],)); rownames(', agList[i], ') <- 1:nrow(', agList[i], ')', sep = '')))
#   return(eval(parse(text = agList[i])))
#   
# })

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

saveRDS(prec,  "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/chirps/prec_new.rds")
saveRDS(srad, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_srad/srad_new.rds")
saveRDS(tmax, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmax/tmax_new.rds")
saveRDS(tmin, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/input_tables/agmerra_tmin/tmin_new.rds")

rm(eli,eli1,prec,ggcmi_data)




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



