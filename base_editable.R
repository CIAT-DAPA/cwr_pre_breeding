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
  continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
  for(k in 1:length(continentList)){
    calc_wat_bal(continent = continentList[k], ncores = 15)
  }; rm(k)
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
    generalIndices(crop = crop, continent = continentList[2], ncores = 20)
  }; rm(k) # Until here
}

# Calculating specific bean indices: loop by continent (improve it)

if(!file.exists()){
  setwd("/home/bmora/cwr_pre_breeding")
  source("02_calculate_index_been_heat.R") # Please run since here
  continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
  for(k in 1:length(continentList)){
    beanIndices(continent = continentList[5], ncores = 15)
  }; rm(k) # Until here
}

# Calculating similarity measures for general crop indices (improve it)
if(!file.exists()){
  source("")
}

# Calculating similarity measures for specific crop indices (improve it)

# Future analysis

###########general index 

crop <- "Finger_millet"
setwd("/home/bmora/cwr_pre_breeding")
source("01_calculating_general_indices.R") # Please run since here
continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
generalIndices(crop = crop, continent = continentList[2], ncores =10)



######### specific index been 

  setwd("/home/bmora/cwr_pre_breeding")
  source("02_calculate_index_been_heat.R") # Please run since here
  continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
     beanIndices(continent = continentList[5], ncores = 10)
  


     ######### specific index Sunflower 
     
     setwd("/home/bmora/cwr_pre_breeding")
     source("02_calculate_index_been_heat.R") # Please run since here
     continentList <- c("Africa", "America", "Asia", "Europa", "Oceania")
     sunflowerIndices(continent = continentList[3], ncores = 10)
     








