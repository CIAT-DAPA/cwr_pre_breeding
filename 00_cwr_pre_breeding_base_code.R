# CWR pre-breeding characterising testing environments: base code
# Authors: B. Mora & H. Achicanoy
# CIAT, 2017

## Bias-correction & downscaling for future data

options(warn = -1); options(scipen = 999); g <- gc(reset = T); rm(list = ls())

# loop by crop
crop <- "bean"

# Step 1: Crop area identification
source("01_cwr_id_occurrence_data.R")
system.time(exp = {occ_data <- crop_area_id(crop = crop); rm(crop_area_id)})

# Step 2: Crop cycle identification by pixel
source("02_cwr_id_crop_growing_cycle.R")
system.time(exp = {occ_data <-crop_cycle_id(crop = "potato", occ_data = occ_data);rm(crop_cycle_id)})

# Step 3: Extract current climate data
source("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Scripts/extractCurrentClimate.R")
current_climate(crop = "potato", occ_data = occ_data)

# Step 4: Index calculation for current data
source("")

# Step 5: Extract future climate
source("")

# Step 6: Index calculation for future data
source("")

# Step 7: Calculate similarity index according to pre-breeding sites
source("")