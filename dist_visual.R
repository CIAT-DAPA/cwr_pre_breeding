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


###############################################################################

tabla2 <-readRDS( paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/global_index_estandarizados_frijol_calor.rds"))
indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada 
cellIdentified <- cellFromXY(base, xy= c(-74.224306,10.727389)) ## Coordenada a buscar
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
saveRDS(tdist,paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/dist_global_index_calor_ciat_palmira.rds"))

##############    FILTRO CULTIVOS   #######################


bean  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Precense_data/area_base.rds"))
#potato  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Potato/database/Precense_data/area_base.rds"))
#finger  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds"))
#eggplant  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Eggplant/database/area_base.rds"))
tdist <- readRDS(paste0(root,"CWR_pre-breeding/Input_data/presence_data/Bean/database/General_indexes/dist_global_index_ciat_palmira.rds"))
crop <- bean
#crop <-  potato 
#crop <-  finger
#crop <-  eggplant


cell <- crop$cellID
tdist2 <- dplyr::filter(tdist, tdist$cellID %in% cell)
dist2 <- tdist2
rm(tdist,tdist2,crop,cell)

#####
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/shape_world"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
shapefile_df <- fortify(countries)
#dist2 <- bean

####################DISTANCIAS DTW
##### Deciles de la distribucion
dist2$DTWarp_cat <- cut(dist2$DTWarp, quantile(dist2$DTWarp, prob = seq(0, 0.5, by = .1)))
dist2$DTWarp_cat <- as.character(dist2$DTWarp_cat)

library(gtools)

qntls <- as.character(na.omit(unique(dist2$DTWarp_cat)))
qntls <- qntls[gtools::mixedorder(qntls)]
pcnts <- paste0(seq(10, 50, 10), "%")

for(i in 1:length(qntls)){
  dist2$DTWarp_cat[which(dist2$DTWarp_cat == qntls[i])] <- pcnts[i]
}; rm(i)

dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
dist2$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist2$DTWarp_cat))


# colours <- c("red", "yellow", "limegreen","royalblue3", "turquoise1")
# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
# colours <- c("red", "yellow", "limegreen","royalblue3", "turquoise1")
# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
marker = list(color = brewer.pal(5, "Set1"))
marker
colours <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")
Y <-ggplot(data = dist2, aes(x = x, y = y, fill = DTWarp_cat)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=shapefile_df, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  coord_cartesian(ylim = c(-50, 50)) +
  ggtitle("Distance DTWarp")+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Similarity (%)"))
Y <- Y + geom_point(aes(x = -76.35666666, y = 3.50472222), shape = 25, colour = "lawngreen", fill = "lawngreen", size = 2)

ggsave(filename=paste0(root, '/CWR_pre-breeding/Input_data/presence_data/Bean/plots/DTW_archivo_calor_palmira.png'),plot=Y, width=10, height=3.5, units='in')


Y <-ggplot(data = dist2, aes(x = x, y = y, fill = DTWarp_cat)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=shapefile_df, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  coord_cartesian(xlim = c(-84, -61), ylim = c(-5, 15)) +
  ggtitle("Distance DTWarp")+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Similarity (%)"))
Y <- Y + geom_point(aes(x = -76.35666666, y = 3.50472222), shape = 25, colour = "lawngreen", fill = "lawngreen", size = 3)
ggsave(filename=paste0(root, '/CWR_pre-breeding/Input_data/presence_data/Bean/plots/DTW_archivo_calor_palmira_zoom.png'),plot=Y, width=8, height=8, units='in')
