# CWR pre-breeding characterising testing environments: master code
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

# Calculating current water balance for all pixels over the world
if(!file.exists()){
  source("00_calculating_current_water_balance.R")
}

# Calculating future water balance for all pixels over the world
if(!file.exists()){
  source("00_calculating_future_water_balance.R")
}

# Calculating general crop indices: loop by continent
crop <- "Bean"
if(!file.exists()){
  source("01_calculating_general_indices.R") # Please run since here
  continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
  for(k in 1:length(continentList)){
    generalIndices(crop = crop, continent = continentList[k], ncores = 15)
  }; rm(i) # Until here
}

# Calculating specific crop indices: loop by continent (improve it)

# Calculating similarity measures for general crop indices (improve it)

# Calculating similarity measures for specific crop indices (improve it)

# Future analysis
