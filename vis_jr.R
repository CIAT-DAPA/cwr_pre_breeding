#CWR pre-breeding characterising testing environments: Index sunflower 
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


rcp <- "rcp85"
gcm <- "gcm1"
crop <- "Bean"
lugar<- "ciat_palmira"



tabla2<- readRDS (paste0(root, "/CWR_pre-breeding/Results/",crop,"/_future/",rcp,"/",gcm,"/Crop_index/global_index_estandarizados_",tolower(crop),"_calor.rds"))
tabla2 <- tabla2[!duplicated(tabla2$cellID),]
indexes <- tabla2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
###coordenada 


CIAT PALMIRA 	 -76,3566666666666	3,50472220000000	
CORPOICA CATIBIA  -74,224306	10,727389
India  78,47416666667	17,375277800
################



cellIdentified <- cellFromXY(base, xy= c(-76.3566666666666 ,	3.50472220000000)) ## Coordenada a buscar
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

saveRDS(tdist,paste0(root, "/CWR_pre-breeding/Results/",crop,"/_future/",rcp,"/",gcm,"/Crop_index/dist_global_index_calor_",lugar,".rds"))
##############    FILTRO CULTIVOS   #######################

tdist <- readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/_future/",rcp,"/",gcm,"/Crop_index/dist_global_index_calor_",lugar,".rds"))
dist2 <- tdist


################# raster con 3 valores
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/_world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
shapefile_df <- fortify(countries)
#dist2 <- bean

####################DISTANCIAS DTW
##### Deciles de la distribucion
dist2$DTWarp_cat <- cut(dist2$DTWarp, quantile(dist2$DTWarp, prob = seq(0,1, by = .1)))
dist2$DTWarp_cat <- as.character(dist2$DTWarp_cat)

library(gtools)
qntls <- as.character(na.omit(unique(dist2$DTWarp_cat)))
qntls <- qntls[gtools::mixedorder(qntls)]
pcnts <- c(paste0(seq(10, 100, 10), "%"))

for(i in 1:length(qntls)){
  dist2$DTWarp_cat[which(dist2$DTWarp_cat == qntls[i])] <- pcnts[i]
}; rm(i)


dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
dist2$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist2$DTWarp_cat))
na <- dist2[which(dist2$DTWarp == 0),]
na.omit(dist2$DTWarp_cat)
na <- data.frame(cellID= na$cellID, DTWarp =1 ,DTWarp_cat= "10%",x= na$x, y= na$y,Cat_num= 10 )
dist2 <- rbind(dist2,na)



# dist2$DTWarp_cat <- as.factor(dist2$DTWarp_cat)
# qntls <- qntls[gtools::mixedorder(qntls)]

dist2$DTWarp_cat <- factor(dist2$DTWarp_cat, levels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))


dist2$DTWarp_cat[which(dist2$DTWarp_cat =="10%")] <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="20%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="30%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="40%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="50%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="60%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="70%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="80%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="90%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="100%")] <- "Low"




colours <- c("#31a354", "#addd8e", "#f7fcb9")
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
Y <- Y + geom_point(aes(dist2[cor,]$x , dist2[cor,]$y ), shape = 25, colour = "#a50f15", fill = "#a50f15", size = 2)

ggsave(filename=paste0(root, "/CWR_pre-breeding/Results/",crop,"/plots/_future/",rcp,"/",gcm,"/dist_global_index_calor_",lugar,".png"),plot=Y, width=10, height=3.5, units='in')

Y <-ggplot(data = dist2, aes(x = x, y = y, fill = DTWarp_cat)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=shapefile_df, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  coord_cartesian(xlim = c(70,90), ylim = c(5, 32)) +
  ggtitle("Distance DTWarp")+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Similarity (%)"))
Y <- Y + geom_point(aes(dist2[cor,]$x , dist2[cor,]$y ), shape = 25, colour = "#a50f15", fill = "#a50f15", size = 5)
ggsave(filename=paste0(root, "/CWR_pre-breeding/Results/",crop,"/plots/_future/",rcp,"/",gcm,"/dist_global_index_calor_",lugar,"_zoom.png"),plot=Y, width=8, height=8, units='in')

#######################################################################################################################################################################
################# UN SOLO MAPA POR STRESS 
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

crop <- "Bean"
lugar<- "corpoica"
tdist2 <- readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/dist_global_index_calor_",lugar,".rds"))
dist2 <- tdist2


indexes <- dist2
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL


# CIAT PALMIRA 	 -76,3566666666666	3,50472220000000	
# CORPOICA CATIBIA  -74,224306	10,727389
# India  78,47416666667	17,375277800

cellIdentified <- cellFromXY(base, xy= c(-76.3566666666666 ,	3.50472220000000)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)

###coordenada 


dist2$DTWarp_cat <- cut(dist2$DTWarp, quantile(dist2$DTWarp, prob = seq(0,1, by = .1)))
dist2$DTWarp_cat <- as.character(dist2$DTWarp_cat)

library(gtools)
qntls <- as.character(na.omit(unique(dist2$DTWarp_cat)))
qntls <- qntls[gtools::mixedorder(qntls)]
pcnts <- c(paste0(seq(10, 100, 10), "%"))

for(i in 1:length(qntls)){
  dist2$DTWarp_cat[which(dist2$DTWarp_cat == qntls[i])] <- pcnts[i]
}; rm(i)


dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
dist2$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist2$DTWarp_cat))


dist2[cor,]$DTWarp <- 1
dist2[cor,]$DTWarp_cat <- "10%"
dist2[cor,]$Cat_num <- 10


dist2$DTWarp_cat <- factor(dist2$DTWarp_cat, levels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))


###############  EL MENOR VALOR 


crop <- "Bean"
lugar<- "ciat_palmira"
tdist3<- readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/dist_global_index_calor_",lugar,".rds"))
dist3 <- tdist3


indexes <- dist3
indexes$Variable <- NULL
require(tidyr)
rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL


CIAT PALMIRA 	 -76.3566666666666 ,	3.50472220000000	
CORPOICA CATIBIA  -74,224306	10,727389
India  78,47416666667	17,375277800

cellIdentified <- cellFromXY(base, xy= c(-76.3566666666666 ,	3.50472220000000)) ## Coordenada a buscar
cor<- which(rownames(indexes) == cellIdentified)

###coordenada 


dist3$DTWarp_cat <- cut(dist3$DTWarp, quantile(dist3$DTWarp, prob = seq(0,1, by = .1)))
dist3$DTWarp_cat <- as.character(dist3$DTWarp_cat)

library(gtools)
qntls <- as.character(na.omit(unique(dist3$DTWarp_cat)))
qntls <- qntls[gtools::mixedorder(qntls)]
pcnts <- c(paste0(seq(10, 100, 10), "%"))

for(i in 1:length(qntls)){
  dist3$DTWarp_cat[which(dist3$DTWarp_cat == qntls[i])] <- pcnts[i]
}; rm(i)


dist3 <- cbind(dist3, xyFromCell(object = base, cell = dist3$cellID))
dist3$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist3$DTWarp_cat))


dist3[cor,]$DTWarp <- 1
dist3[cor,]$DTWarp_cat <- "10%"
dist3[cor,]$Cat_num <- 10


dist3$DTWarp_cat <- factor(dist3$DTWarp_cat, levels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))

####### CALCULO DEL MINIMO 
dist2$Cat_num1 <- dist3$Cat_num
tabla <- dist2
tabla$Cat_numr<- dat <- transform(tabla, min = pmin(Cat_num,Cat_num1))
min <- tabla[,ncol(tabla)]
dist2 <- min 
dist2$Cat_num1 <- NULL
dist2$Cat_num <- NULL
dist2$min <- dist2$min /100


library(formattable)
dist2$DTWarp_cat<- as.character(percent(dist2$min))
dist2$min <- NULL
dist2$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist3$DTWarp_cat))

d <- dist2
#
dist2<- d

dist2$DTWarp_cat <- factor(dist2$DTWarp_cat, levels = c("10.00%","20.00%","30.00%","40.00%","50.00%","60.00%","70.00%","80.00%","90.00%","100.00%"))
dist2$DTWarp_cat <- as.character(dist2$DTWarp_cat)



dist2[cor,]$DTWarp <- 1
dist2[cor,]$DTWarp_cat <- "10.00%"
dist2[cor,]$Cat_num <- 10





dist2$DTWarp_cat[which(dist2$DTWarp_cat =="10.00%")] <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="20.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="30.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="40.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="50.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="60.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="70.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="80.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="90.00%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="100.00%")] <- "Low"

na.omit(dist2$DTWarp_cat)


################# raster con 3 valores
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/_world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
countries$COUNTRY <- iconv(countries$COUNTRY, from = "UTF-8", to = "latin1")
shapefile_df <- fortify(countries)
#dist2 <- bean


colours <- c("#31a354", "#addd8e", "#f7fcb9")
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
Y <- Y + geom_point(aes(dist2[cor,]$x , dist2[cor,]$y ), shape = 25, colour = "#a50f15", fill = "#a50f15", size = 2)

ggsave(filename=paste0(root, "/CWR_pre-breeding/Results/",crop,"/plots/_future/",rcp,"/",gcm,"/dist_global_index_calor_",lugar,".png"),plot=Y, width=10, height=3.5, units='in')

Y <-ggplot(data = dist2, aes(x = x, y = y, fill = DTWarp_cat)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=shapefile_df, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  coord_cartesian(xlim = c(70,90), ylim = c(5, 32)) +
  ggtitle("Distance DTWarp")+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Similarity (%)"))
Y <- Y + geom_point(aes(dist2[cor,]$x , dist2[cor,]$y ), shape = 25, colour = "#a50f15", fill = "#a50f15", size = 5)
ggsave(filename=paste0(root, "/CWR_pre-breeding/Results/",crop,"/plots/_future/",rcp,"/",gcm,"/dist_global_index_calor_",lugar,"_zoom.png"),plot=Y, width=8, height=8, units='in')

