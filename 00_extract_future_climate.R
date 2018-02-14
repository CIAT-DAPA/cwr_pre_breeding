library(raster)
library(sp)
library(rgdal)
library(velox)
library(tidyverse)

avlb_rst <- list.files(path = "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k", pattern = "*.nc$", full.names = F)
avlb_rst <- avlb_rst[setdiff(1:length(avlb_rst), grep(pattern = "historical", x = avlb_rst))]
avlb_rst <- avlb_rst[-1]

test <- data.frame(do.call('rbind', strsplit(as.character(avlb_rst),'_',fixed=TRUE)))
colnames(test) <- c("Variable", "Time", "D1", "RCP", "D2", "GCM", "Year")

avlb_rst[setdiff(1:length(avlb_rst), grep(pattern = "^pr_day_", x = avlb_rst))]

# Extract data from one year
tmax <- raster::stack("//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k/tasmax_day_BCSD_rcp45_r1i1p1_MRI-CGCM3_2043.nc")
tmax_vx <- velox::velox(tmax)
pnt <- tmax_vx$getCoordinates()
colnames(pnt) <- c("x", "y")
pnt <- data.frame(pnt)
pnt$x[pnt$x > 180] <- pnt$x[pnt$x > 180] - 360
rtmp <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/AgMerra_template.RDS")
rtmp[which(!is.na(rtmp[]))] <- 1
# pnt_rtmp <- coordinates(rtmp[which(!is.na(rtmp[]))])
pnt_rtmp <- rasterToPoints(rtmp)
pnt_rtmp <- pnt_rtmp[,1:2]
pnt_rtmp <- data.frame(pnt_rtmp)
pnt_rtmp$cellID <- cellFromXY(object = rtmp, xy = pnt_rtmp)
pnt$cellID <- cellFromXY(object = rtmp, xy = pnt)
pnt <- pnt[pnt$cellID %in% pnt_rtmp$cellID,]
rownames(pnt) <- 1:nrow(pnt)
pnt$x[pnt$x <= 0] <- pnt$x[pnt$x <= 0] + 360
vls <- tmax_vx$extract_points(sp = SpatialPoints(pnt[,1:2]))
vls <- as_data_frame(vls)
vls <- cbind(pnt, vls)
vls$x[vls$x > 180] <- vls$x[vls$x > 180] - 360

