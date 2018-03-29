# Brayan Mora


############################################################################################################################################
############################################################################################################################################
##############################################      85                   ####################################################################
############################################################################################################################################
############################################################################################################################################

future <- list.files(path = "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k", pattern = "*.nc$", full.names = F)


#gcm1_prec <-  future[grep(pattern = "tasmax_day_BCSD_rcp85_r1i1p1_bcc-csm1-1", x = future)]
# gcm1_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_bcc-csm1-1", x = future)]
# gcm1_prec <-  future[grep(pattern = "pr_day_BCSD_rcp85_r1i1p1_bcc-csm1-1", x = future)]


###########################################################################################################
#gcm2_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_MIROC-ESM-CHEM", x = future)]
#gcm2_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_MIROC-ESM-CHEM", x = future)]
gcm2_prec <-  future[grep(pattern = "tasmax_day_BCSD_rcp85_r1i1p1_MIROC-ESM-CHEM", x = future)]




###########################################################################################################
gcm3_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_BNU-ESM", x = future)]
gcm3_prec <-  future[grep(pattern = "tasmax_day_BCSD_rcp85_r1i1p1_BNU-ESM", x = future)]
gcm3_prec <-  future[grep(pattern = "pr_day_BCSD_rcp85_r1i1p1_BNU-ESM", x = future)]



###########################################################################################################
gcm4_prec <-  future[grep(pattern = "tasmax_day_BCSD_rcp85_r1i1p1_MPI-ESM-LR", x = future)]
gcm4_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_MPI-ESM-LR", x = future)]
gcm4_prec <-  future[grep(pattern = "pr_day_BCSD_rcp85_r1i1p1_MPI-ESM-LR", x = future)]

###########################################################################################################
#gcm5_prec <-  future[grep(pattern = "tasmin_day_BCSD_rcp85_r1i1p1_inmcm4.", x = future)]
gcm5_prec <-  future[grep(pattern = "tasmax_day_BCSD_rcp85_r1i1p1_inmcm4.", x = future)]


g <- gcm5_prec[-(1:14)]
g1 <- g[-(51:length(g))]


library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

cores<- detectCores()
cl<- makeCluster(cores-20)
registerDoParallel(cl) 


# mylist <- rep(list(NA),50)

system.time(l <- foreach(i=1:50) %dopar% {
  mylist <- list()
  root <- "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k/"
  require(dplyr)
  require(velox)
  require(raster)
  require(raster)
  
  cat(paste0("procesando tabla:",i,"\n" ))
  prec <- raster::stack(paste0(root, g[i]))  
  prec_vx <- velox::velox(prec)
  pnt <- prec_vx$getCoordinates()
  colnames(pnt) <- c("x", "y")
  pnt <- data.frame(pnt)
  pnt$x[pnt$x > 180] <- pnt$x[pnt$x > 180] - 360
  rtmp <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/AgMerra_template.RDS")
  rtmp[which(!is.na(rtmp[]))] <- 1
  # rtmp[which(!is.na(rtmp[]))] <- 1
  # pnt_rtmp <- coordinates(rtmp[which(!is.na(rtmp[]))])
  pnt_rtmp <- rasterToPoints(rtmp)
  pnt_rtmp <- pnt_rtmp[,1:2]
  pnt_rtmp <- data.frame(pnt_rtmp)
  pnt_rtmp$cellID <- cellFromXY(object = rtmp, xy = pnt_rtmp)
  pnt$cellID <- cellFromXY(object = rtmp, xy = pnt)
  pnt <- pnt[pnt$cellID %in% pnt_rtmp$cellID,]
  rownames(pnt) <- 1:nrow(pnt)
  pnt$x[pnt$x <= 0] <- pnt$x[pnt$x <= 0] + 360
  vls <- prec_vx$extract_points(sp = SpatialPoints(pnt[,1:2]))
  vls <- as_data_frame(vls)
  vls <- cbind(pnt, vls)
  vls$x[vls$x > 180] <- vls$x[vls$x > 180] - 360
  vls1<- data.frame(cellID= vls$cellID, lon= vls$x, lat= vls$y, vls[,-3])
  vls1$x <- NULL 
  vls1$y <- NULL
  names(vls1)[4:ncol(vls1)] <- names(prec)
  names(vls1)[4:ncol(vls1)] <- as.character(gsub(pattern = ".", replacement = "-", x = gsub(pattern = "X", replacement = "", x = names(vls1)[4:ncol(vls1)]), fixed = T))
  mylist[[i]] <- vls1
  
  
} )
stopCluster(cl)

l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:3)]
  return(x)
})
tmin  <- do.call(cbind,l1)
tmin <- cbind(l[[1]], tmin)


##Filtro

america <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_america.rds")
africa  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_africa.rds")
oceania <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_oceania.rds")
asia    <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_asia.rds")
europa  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_current_climate/agmerra_tmin/tmin_filtered_europa.rds")



a <- list(america,africa,oceania,asia,europa)

filtro <- lapply(1:length(a), function(i){
  tabla <- dplyr::filter(tmin, cellID %in% a[[i]]$cellID)
  tabla1 <- tabla[,-(1:3)]
  t <- tabla1 -273
  tab <- cbind(tabla[,(1:3)], t)
  return(tab)
})

dapa <-"//dapadfs/Workspace_cluster_9"


saveRDS(filtro[[1]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm5/agmerra_tmax/tmax_filtered_america.rds"))
saveRDS(filtro[[2]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm5/agmerra_tmax/tmax_filtered_africa.rds"))
saveRDS(filtro[[3]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm5/agmerra_tmax/tmax_filtered_oceania.rds"))
saveRDS(filtro[[4]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm5/agmerra_tmax/tmax_filtered_asia.rds"))
saveRDS(filtro[[5]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm5/agmerra_tmax/tmax_filtered_europa.rds"))

##############prec


# saveRDS(filtro[[1]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm4/chirps/prec_filtered_america.rds"))
# saveRDS(filtro[[2]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm4/chirps/prec_filtered_africa.rds"))
# saveRDS(filtro[[3]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm4/chirps/prec_filtered_oceania.rds"))
# saveRDS(filtro[[4]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm4/chirps/prec_filtered_asia.rds"))
# saveRDS(filtro[[5]], paste0(dapa,"/CWR_pre-breeding/Input_data/_future_climate/rcp85/gcm4/chirps/prec_filtered_europa.rds"))
# # # 
# # # 
