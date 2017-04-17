# CWR pre-breeding characterising testing environments: base code
# Authors: B. Mora & H. Achicanoy
# CIAT, 2017

## Bias-correction & downscaling for future data

options(warn = -1); options(scipen = 999); g <- gc(reset = T); rm(list = ls())

suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})

# loop by crop

# Step 1: Crop area identification
source("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Scripts/Datos_presencia_potato.R")
crop_area_id(crop = "potato")

# Step 2: Crop cycle identification by pixel
source("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Scripts/Ciclo_cultivo_potato.R")
crop_cycle_id(crop = "potato", sys_type = "rainfed") # rainfed, irrigated

# Step 3: Extract current climate data
source("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Scripts/extractCurrentClimate.R")
current_climate(crop = "potato")

# Step 4: Index calculation for current data
source("")

# Step 5: Extract future climate
source("")

# Step 6: Index calculation for future data
source("")

# Step 7: Calculate similarity index according to pre-breeding sites
source("")