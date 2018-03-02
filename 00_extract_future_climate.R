# CWR pre-breeding characterising testing environments: extract future climate
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

# Load R packages
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(sp)){install.packages('sp'); library(sp)} else {library(sp)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(velox)){install.packages('velox'); library(velox)} else {library(velox)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})

# # Matching combinations GCM-RCP-time
# c("BNU-ESM" = "bnu_esm")
# c("CanESM2" = "cccma_canesm2")
# c("CSIRO-Mk3-6-0" = "csiro_mk3_6_0")
# c("GFDL-CM3" = "gfdl_cm3")
# c("inmcm4" = "inm_cm4")
# c("IPSL-CM5A-LR" = "ipsl_cm5a_lr")
# c("IPSL-CM5A-MR" = "ipsl_cm5a_mr")
# c("MIROC-ESM-CHEM" = "miroc_esm_chem")
# c("MPI-ESM-LR" = "mpi_esm_lr")
# c("MPI-ESM-MR" = "mpi_esm_mr")
# c("NorESM1-M" = "ncc_noresm1_m")

climateWizardSrc <- "//ccafsserver.cgiarad.org/data_climatewizard/AR5_Global_Daily_25k" # It has: prec, tmax, tmin
rstList <- list.files(path = climateWizardSrc, pattern = "*.nc$", full.names = F)
rstList <- rstList[setdiff(1:length(rstList), grep(pattern = "historical", x = rstList))]
rstList <- rstList[-1]
grep2 <- Vectorize(FUN = grep, vectorize.args = "pattern")
rstList <- rstList[grep2(pattern = paste(2020:2069), x = rstList)]
rstList <- rstList[grep2(pattern = c("BNU-ESM", "CanESM2", "CSIRO-Mk3-6-0",
                                     "GFDL-CM3", "inmcm4", "IPSL-CM5A-LR",
                                     "IPSL-CM5A-MR", "MIROC-ESM-CHEM", "MPI-ESM-LR",
                                     "MPI-ESM-MR", "NorESM1-M"), x = rstList)]
df <- data.frame(do.call('rbind', strsplit(as.character(rstList), '_', fixed = TRUE)))
colnames(df) <- c("Variable", "Time", "D1", "RCP", "D2", "GCM", "Year")
gcmList <- df$GCM %>% as.character %>% unique

ownSrc <- "//dapadfs/workspace_cluster_12/Trust"
own_rcp45_gcmList <- list.dirs(path = paste0(ownSrc, "/gcm_rcp45/bc_0_5deg_lat"), full.names = F, recursive = F)
own_rcp85_gcmList <- list.dirs(path = paste0(ownSrc, "/gcm_rcp85/bc_0_5deg_lat"), full.names = F, recursive = F)





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

