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

america <-  readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/index_tabla_rbind_frijol_america.rds"))
africa<-    readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/index_tabla_rbind_frijol_africa.rds"))
asia<-      readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/index_tabla_rbind_frijol_asia.rds"))
europa<-    readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/index_tabla_rbind_frijol_europa.rds"))
oceania<-   readRDS(paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/index_tabla_rbind_frijol_oceania.rds"))

unique(america$Variable)


america<- dplyr::filter(america, Variable == "Texteme" | Variable == "days_tmin22" |Variable == "days_tmin24" |Variable == "calor22"|Variable == "calor24"|Variable == "optimo" )
africa<- dplyr::filter(africa, Variable == "Texteme" | Variable == "days_tmin22" |Variable == "days_tmin24" |Variable == "calor22"|Variable == "calor24"|Variable == "optimo" )
europa<- dplyr::filter(europa, Variable == "Texteme" | Variable == "days_tmin22" |Variable == "days_tmin24" |Variable == "calor22"|Variable == "calor24"|Variable == "optimo" )
asia<- dplyr::filter(asia,Variable == "Texteme" | Variable == "days_tmin22" |Variable == "days_tmin24" |Variable == "calor22"|Variable == "calor24"|Variable == "optimo" )
oceania<- dplyr::filter(oceania,Variable == "Texteme" | Variable == "days_tmin22" |Variable == "days_tmin24" |Variable == "calor22"|Variable == "calor24"|Variable == "optimo" )


#tabla <-  america
#tabla <-  africa
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
#africa <- tabla2
#europa <- tabla2
##asia <- tabla2
oceania <- tabla2
rm(tabla,tabla2)

tabla2 <- rbind(america, europa,asia,africa,oceania)
saveRDS(tabla2, paste0(root, "/CWR_pre-breeding/Input_data/presence_data/Bean/database/Crop_indexes/global_index_estandarizados_frijol_calor.rds"))
rm(america,africa,asia,oceania,europa)




