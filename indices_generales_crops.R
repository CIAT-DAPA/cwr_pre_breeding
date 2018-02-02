# Load packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
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
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
suppressMessages(if(!require(googlesheets)){install.packages("googlesheets");library(googlesheets)}else{library(googlesheets)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})

# Path options
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

cat("Verify if Standarized global indices exist\n")
if(!file.exists(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/global_index_estandarizados.rds"))){
  
  cat("Standarized global indices exist don't exist. Calculating ...\n")
  # Load general indices by continent
  america <- readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_america_completo.rds"))
  africa <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_africa_completo.rds"))
  europa <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_europa_completo.rds"))
  asia <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_asia_completo.rds"))
  oceania <- readRDS(paste0(root,"/CWR_pre-breeding/Results/input_tables/chirps/index_tabla_rbind_oceania_completo.rds"))
  
  # Select specific indices
  america<- dplyr::filter(america, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
  africa<- dplyr::filter(africa, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
  europa<- dplyr::filter(europa, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
  asia<- dplyr::filter(asia, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
  oceania<- dplyr::filter(oceania, Variable == "TOTRAIN" | Variable == "CDD" |Variable == "P5D" |Variable == "P_95" )
  
  continents <- c("america", "africa", "europa", "asia", "oceania")
  process_continents <- lapply(1:length(continents), function(i){
    
    eval(parse(text = paste('tabla <- ', continents[i], sep = '')))
    
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
    
    return(tabla2)
    
  })
  
  tabla2 <- do.call(rbind, process_continents)
  saveRDS(tabla2, paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/global_index_estandarizados.rds"))
  rm(america, africa, asia, oceania, europa, process_continents)
  
} else {
  
  cat("Calculate similarity index by crop\n")
  
  tabla2 <-readRDS(paste0(root, "/CWR_pre-breeding/Results/input_tables/chirps/global_index_estandarizados.rds"))
  indexes <- tabla2
  indexes$Variable <- NULL
  rownames(indexes) <- indexes$cellID; indexes$cellID <- NULL
  
  # Load database from Google Drive
  crop_info <- gs_ls("crop_information_states")
  crop_info <- gs_title("crop_information_states")
  crop_info %>% gs_browse(ws = "Status")
  crop_info <- crop_info %>% gs_read(ws = "Status")
  crop_info$Longitude <- gsub(pattern = ",", replacement = ".", x = crop_info$Longitude) %>% as.character %>% as.numeric
  crop_info$Latitude <- gsub(pattern = ",", replacement = ".", x = crop_info$Latitude) %>% as.character %>% as.numeric
  
  cropList <- unique(crop_info$Crop) %>% as.character
  
  crop_site_similarity <- function(crop = "Bean"){
    
    crop_data <- crop_info %>% filter(Crop == crop)
    
    cat("Identifying coordinates and calculating similarity indices\n")
    lapply(1:nrow(crop_data), function(i){
      
      if(!file.exists(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/", crop, "/database/dist_global_index_", crop_data$Place2[i], ".rds"))){
        
        cat("Selecting coordinate:", i,"\n")
        cellIdentified <- cellFromXY(base, xy = c(crop_data$Longitude[i], crop_data$Latitude[i]))
        cor <- which(rownames(indexes) == cellIdentified)
        tabla2[is.na(tabla2)] <- 0
        
        cat("Calculating DTW similarity index\n")
        dissMat <- matrix(data = NA, nrow = nrow(tabla2), ncol = 1)
        system.time (for(i in 1:nrow(tabla2)){
          dissMat[i,1] <- TSclust::diss.DTWARP(x = as.numeric(tabla2[i, -1]), y = as.numeric(tabla2[cor, -1]))
        }); rm(i)
        
        dissMat <- as.data.frame(dissMat)
        colnames(dissMat) <- c("DTWarp")
        tdist <- data.frame(cellID = tabla2$cellID, dissMat)
        rm(dissMat)
        saveRDS(tdist, paste0(root, "/CWR_pre-breeding/Input_data/presence_data/", crop, "/database/dist_global_index_", crop_data$Place2[i], ".rds"))
        
      } else {
        cat("Similarity index for:", crop_data$Place[i], " exists. Don't require to calculate again.\n")
      }
      
    })
    
    return(cat("Process done successfully\n"))
  }
  lapply(X = cropList, function(x) crop_site_similarity(crop = x))
}

########################################################################################################
########################################################################## 
##############        FILTRO CULTIVOS    #################################
########################################################################## 
#####


bean  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds"))
#potato  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Potato/database/area_base.rds"))
#finger  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/area_base.rds"))
#eggplant  <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Eggplant/database/area_base.rds"))

tdist <- readRDS(paste0(root,"/CWR_pre-breeding/Input_data/presence_data/Bean/database/dist_global_index_palmira.rds"))

crop <- bean
#crop <-  potato 
#crop <-  finger
#crop <-  eggplant


cell <- crop$cellID
tdist2 <- dplyr::filter(tdist, tdist$cellID %in% cell)
dist2 <- tdist2
rm(tdist,tdist2,crop,cell)

#####
countries <- rgdal::readOGR(dsn = paste0(root, "/CWR_pre-breeding/Input_data/world_shape"), "all_countries")
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
rm(tdist2)

base[] <- NA
base[dist2$cellID] <- dist2$Cat_num
plot(base)
writeRaster(x = base, filename = paste0(root, "/CWR_pre-breeding/Results/Deciles_bean.tif"))

library(leaflet)
crs(base) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(base), na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(base, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(base),
            title = "Similarity locations to CIAT-Palmira")



# colours <- c("red", "yellow", "limegreen","royalblue3", "turquoise1")
# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
marker = list(color = brewer.pal(5, "Set1"))
marker
colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
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
Y <- Y + geom_point(aes(x = -76.356666, y = 3.5047222), shape = 25, colour = "lawngreen", fill = "lawngreen", size = 2)

ggsave(filename="DTWarp_archivo2.png", plot=Y, width=10, height=3.5, units='in')

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
Y <- Y + geom_point(aes(x = -76.356666, y = 3.5047222), shape = 25, colour = "lawngreen", fill = "lawngreen", size = 2)
ggsave(filename="DTWarp_archivo3.png", plot=Y, width = 8, height = 8, units='in')

ggsave(filename=paste0(root, '/CWR_pre-breeding/Input_data/presence_data/Eggplant/plots/DTW_archivo_valencia.png'), plot=Y, width=10, height=8, units='in')





