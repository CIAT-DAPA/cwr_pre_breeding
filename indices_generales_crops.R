
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

america <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_america_completo.rds"))
africa <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_africa_completo.rds"))
europa<- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_europa_completo.rds"))
asia<- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_asia_completo.rds"))
oceania <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_oceania_completo.rds"))

america<- dplyr::filter(america, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
africa<- dplyr::filter(africa, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
europa<- dplyr::filter(europa, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
asia<- dplyr::filter(asia, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
oceania<- dplyr::filter(oceania, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )

# tabla <-  america
# tabla <-  africa
# tabla <-  europa
# tabla <-  asia
tabla <- oceania


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
# america <- tabla2
# africa <- tabla2
# europa <- tabla2
# asia <- tabla2
oceania <- tabla2
rm(tabla,tabla2)

tabla2 <- rbind(america, europa,asia,africa,oceania)
saveRDS(tabla2, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/global_index_estandarizados.rds"))
rm(america,africa,asia,oceania,europa)

########################################################################################################################################
########################################################################################################################################
####################################       LO DE ATRAS NO SE VUELVE A CORRER        ####################################################
########################################################################################################################################
########################################################################################################################################

#cORRER DE AQUI EN ADELANTE PARA CADA CULTIVO 

tabla2 <-readRDS( paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/global_index_estandarizados.rds"))
indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada 
cellIdentified <- cellFromXY(base, xy= c(-0.33416666666666667,39.4777778)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)
tabla2[is.na(tabla2)] <- 0

 
dissMat <- matrix(data = NA, nrow = nrow(tabla2), ncol = 1)
system.time (for(i in 1:nrow(tabla2)){
  dissMat[i,1] <- TSclust::diss.DTWARP(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
}); rm(i)


dissMat <- as.data.frame(dissMat)
colnames(dissMat) <- c( "DTWarp")
require(tidyr)
tdist <- data.frame(cellID = tabla2$cellID, dissMat)
rm(dissMat)
saveRDS(tdist,paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Eggplant/database/dist_global_index_valencia.rds"))

########################################################################################################
########################################################################## 
##############        FILTRO CULTIVOS    #################################
########################################################################## 
#####


# bean  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds"))
# potato  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Potato/database/area_base.rds"))
#finger  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds"))
eggplant  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Eggplant/database/area_base.rds"))

#crop <- bean
#crop <-  potato 
#crop <-  finger
crop <-  eggplant


cell <- crop$cellID
tdist2 <- dplyr::filter(tdist, tdist$cellID %in%cell)
rm(tdist,crop,cell)

#####
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
shapefile_df <- fortify(countries)
dist2 <- tdist2

####################DISTANCIAS DTW
##### Deciles de la distribucion
dist2$DTWarp <- cut(dist2$DTWarp, quantile(dist2$DTWarp, prob = seq(0, 0.5, by = .1)))
dist2$DTWarp <- as.character(dist2$DTWarp)
unique(dist2$DTWarp)
##### El valor de los deciles cambia de acuerdo a la linea anterior
dist2$DTWarp[which(dist2$DTWarp == "(0,92.6]")] <- "10%"
dist2$DTWarp[which(dist2$DTWarp == "(92.6,95]")] <- "20%"
dist2$DTWarp[which(dist2$DTWarp == "(95,96.7]")] <- "30%"
dist2$DTWarp[which(dist2$DTWarp == "(96.7,98]")] <- "40%"
dist2$DTWarp[which(dist2$DTWarp == "(98,99.2]")] <- "50%"
dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
rm(tdist2)


colours <- c("red", "yellow", "limegreen","royalblue3", "turquoise1")
Y <-ggplot(data = dist2, aes(x = x, y = y, fill = DTWarp)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=shapefile_df, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  scale_y_continuous(limits = c(-50, 50))+
  ggtitle("Distance  DTWarp ")+
  theme_bw()

ggsave(filename="DTWarp_archivo.png", plot=Y, width=10, height=8, units='in')
ggsave(filename=paste0(root, '/CWR_pre-breeding/Input_data/presence_data/Eggplant/plots/DTW_archivo_valencia.png'), plot=Y, width=10, height=8, units='in')





