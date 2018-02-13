# CWR pre-breeding characterising testing environments: base code
# Authors: B. Mora & H. Achicanoy
# CIAT, 2018

# R options
options(warn = -1); options(scipen = 999); g <- gc(reset = T); rm(list = ls())

# Extract current climate data
if(!file.exists()){
  source("00_extract_current_climate.R")
}

# Extract future climate data
if(!file.exists()){
  source("00_extract_future_climate.R")
}


# loop by crop
crop <- "bean"

# Step 1: Crop area identification
source("01_cwr_id_occurrence_data.R")
system.time(exp = {occ_data <- crop_area_id(crop = crop); rm(crop_area_id)})

# Step 2: Crop cycle identification by pixel
source("02_cwr_id_crop_growing_cycle.R")
system.time(exp = {occ_data <- crop_cycle_id(crop = crop, occ_data = occ_data); rm(crop_cycle_id)})

# Step 3: Index calculation for current data
source("03_cwr_index_calculation_current_climate.R")

# Step 5: Extract future climate
source("")

# Step 6: Index calculation for future data
source("")

# Step 7: Calculate similarity index according to pre-breeding sites
source("")