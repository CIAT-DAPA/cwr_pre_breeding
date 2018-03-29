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


crop <- "Sunflower" 

africa <-  readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/",tolower(crop),"_crop_indices_africa.rds"))
america<-    readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/",tolower(crop),"_crop_indices_america.rds"))
asia<-      readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/",tolower(crop),"_crop_indices_asia.rds"))
europa<-    readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/",tolower(crop),"_crop_indices_europa.rds"))
oceania<-   readRDS(paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/",tolower(crop),"_crop_indices_asia.rds"))

variable <- c(as.character(unique(america$Variable)))

america<- dplyr::filter(america, Variable == variable[1]| Variable == variable[2] |Variable == variable[3] |Variable == variable[4]|Variable == variable[5]|Variable ==variable[6])
africa<-   dplyr::filter(africa, Variable == variable[1] | Variable == variable[2] |Variable == variable[3] |Variable == variable[4]|Variable == variable[5]|Variable == variable[6] )
europa<-   dplyr::filter(europa, Variable == variable[1] | Variable == variable[2] |Variable == variable[3] |Variable == variable[4]|Variable == variable[5]|Variable == variable[6] )
asia<-   dplyr::filter(asia,Variable == variable[1] |Variable == variable[2] |Variable ==  variable[3] |Variable == variable[4]|Variable ==variable[5]|Variable ==variable[6] )
oceania<-   dplyr::filter(oceania,Variable ==variable[1] | Variable == variable[2]|Variable == variable[3] |Variable == variable[4]|Variable == variable[5]|Variable == variable[6] )


tabla <- list(america, africa, europa, asia,oceania) 
tabla2 <- lapply(1:length(tabla), function(j){
  cat(paste0("Processed tabla:", j, "\n"))
  
  tabla[[j]]$Year <- as.numeric(tabla[[j]]$Year)
  tabla[[j]]$Value <- as.numeric(tabla[[j]]$Value)
  tabla[[j]]$cellID <- as.numeric(tabla[[j]]$cellID)
  vars <- c(unique(tabla[[j]]$Variable))
  
  spread_tables <- lapply(1:length(vars), function(i){
    df <- tabla[[j]] %>% filter(Variable == vars[i])
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
  tabla2 <- replace(tabla2, is.na(tabla2),0)
  return(tabla2)
})

tabla2 <- do.call(rbind, tabla2)
saveRDS(tabla2, paste0(root, "/CWR_pre-breeding/Results/",crop,"/Crop_indices/global_index_estandarizados_",tolower(crop),"_calor.rds"))
rm(america,africa,asia,oceania,europa)
