# Calculating summary metrics by crop-stress
# H. Achicanoy
# CIAT, 2018

# R options
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))

# Path settings
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

# ======================================================================= #
# Functions
# ======================================================================= #

# Function to calculate the ensemble median through all available GCMs
# and the minimum distance when there are several pre-breeding sites
min_ftr_dist <- function(all_gcm = T, f_path = f_path, gcm_path = gcm_path){
  
  # All GCMs
  # path <- "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/Bean/_future/rcp85"
  # One GCM
  # path <- "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/Bean/_future/rcp85/gcm1/Crop_index"
  
  if(all_gcm){
    
    # Identify all distance files
    files <- paste0(f_path, "/gcm", 1:5, "/Crop_index") %>% purrr::map(.f = function(x) list.files(x, pattern = "^dist_global*", full.names = T))
    files <- paste0(f_path, "/gcm", 1:5, "/Crop_index") %>% purrr::map(.f = function(x) list.files(x, pattern = "Carrot_index_drought_europa.rds", full.names = T))
    
    files <- files %>% purrr::discard(~ length(.x) == 0)
    
    # Determine number of coordinates
    ncoor <- length(files[[1]])
    # Calculate median ensemble of distances through all GCMs
    df <- lapply(X = 1:ncoor, function(coor){
      
      coor_info <- lapply(X = 1:5, function(gcm){df <- readRDS(files[[gcm]][coor])})
      coor_info <- Reduce(function(x, y) merge(x, y, by = "cellID"), coor_info)
      colnames(coor_info)[2:ncol(coor_info)] <- paste0("gcm", 1:(ncol(coor_info)-1))
      coor_info$DTWarp <- apply(coor_info[,-1], 1, median)
      coor_info <- coor_info %>% dplyr::select(cellID, DTWarp)
      
      coor_info$DTWarp_cat <- cut(coor_info$DTWarp, quantile(coor_info$DTWarp, prob = seq(0, 1, by = .1)))
      coor_info$DTWarp_cat <- as.character(coor_info$DTWarp_cat)
      qntls <- as.character(na.omit(unique(coor_info$DTWarp_cat)))
      qntls <- qntls[gtools::mixedorder(qntls)]
      pcnts <- paste0(seq(10, 100, 10), "%")
      for(i in 1:length(qntls)){
        coor_info$DTWarp_cat[which(coor_info$DTWarp_cat == qntls[i])] <- pcnts[i]
      }; rm(i)
      coor_info$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", coor_info$DTWarp_cat))
      coor_info$DTWarp_cat <- NULL
      coor_info$DTWarp <- NULL
      coor_info <- coor_info %>% dplyr::select(cellID, Cat_num)
      
      return(coor_info)
      
    })
    if(length(df) > 1){
      df <- do.call("merge", append(df, list(by = "cellID")))
      df$Cat_num <- apply(df[,-1], 1, min)
      df <- df %>% dplyr::select(cellID, Cat_num)
    } else {
      if(length(df) == 1){
        df <- as.data.frame(df[[1]])
        df$Cat_num <- df[,2]
        df <- df %>% dplyr::select(cellID, Cat_num)
      }
    }
    
    return(df)
    
  } else {
    
    files <- list.files(path = gcm_path, pattern = "^dist_global*", full.names = T)
    
    if(length(files) > 1){
      
      df <- files %>% purrr::map(.f = readRDS)
      df <- do.call("merge", append(df, list(by = "cellID")))
      df$DTWarp <- apply(df[,-1], 1, median)
      df <- df %>% dplyr::select(cellID, DTWarp)
      df$DTWarp_cat <- cut(df$DTWarp, quantile(df$DTWarp, prob = seq(0, 1, by = .1)))
      df$DTWarp_cat <- as.character(df$DTWarp_cat)
      qntls <- as.character(na.omit(unique(df$DTWarp_cat)))
      qntls <- qntls[gtools::mixedorder(qntls)]
      pcnts <- paste0(seq(10, 100, 10), "%")
      for(i in 1:length(qntls)){
        df$DTWarp_cat[which(df$DTWarp_cat == qntls[i])] <- pcnts[i]
      }; rm(i)
      df$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", df$DTWarp_cat))
      df$DTWarp_cat <- NULL
      df$DTWarp <- NULL
      df <- df %>% dplyr::select(cellID, Cat_num)
      return(df)
      
    } else {
      
      if(length(files) == 1){
        
        df <- readRDS(files)
        df$DTWarp_cat <- cut(df$DTWarp, quantile(df$DTWarp, prob = seq(0, 1, by = .1)))
        df$DTWarp_cat <- as.character(df$DTWarp_cat)
        qntls <- as.character(na.omit(unique(df$DTWarp_cat)))
        qntls <- qntls[gtools::mixedorder(qntls)]
        pcnts <- paste0(seq(10, 100, 10), "%")
        for(i in 1:length(qntls)){
          df$DTWarp_cat[which(df$DTWarp_cat == qntls[i])] <- pcnts[i]
        }; rm(i)
        df$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", df$DTWarp_cat))
        df$DTWarp_cat <- NULL
        return(df)
        
      } else {
        cat("There is no distance data to process. Please check.\n")
      }
      
    }
    
  }
  
}

# Function to transform distances table to raster
tbl_to_rst <- function(path = path, file = T, df = df, current = T){
  
  if(file){
    df <- df
  } else {
    df <- readRDS(path)
  }
  
  if("Cat_num" %in% colnames(df)){
    df <- df %>% dplyr::filter(Cat_num == 10)
  } else {
    df$DTWarp_cat <- cut(df$DTWarp, quantile(df$DTWarp, prob = seq(0, 1, by = .1)))
    df$DTWarp_cat <- as.character(df$DTWarp_cat)
    
    qntls <- as.character(na.omit(unique(df$DTWarp_cat)))
    qntls <- qntls[gtools::mixedorder(qntls)]
    pcnts <- paste0(seq(10, 100, 10), "%")
    for(i in 1:length(qntls)){
      df$DTWarp_cat[which(df$DTWarp_cat == qntls[i])] <- pcnts[i]
    }; rm(i)
    df$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", df$DTWarp_cat))
    df$DTWarp_cat <- NULL
  }
  
  if(current){
    tmpl <- base
    tmpl[] <- NA
    tmpl[df$cellID] <- 1
  } else {
    tmpl <- base
    tmpl[] <- NA
    tmpl[df$cellID] <- 2
  }
  
  return(tmpl)
  
}

# Plotting potential areas raster
plot_maps <- function(lyr = chg_map, shp_path = shp_path, outfile = outfile){
  
  # lyr = chg_map
  shp <- raster::shapefile(shp_path)
  
  gg2 <- rasterVis::gplot(lyr) +
    geom_tile(aes(fill = as.factor(value))) +
    # scale_colour_hue(na.value = 'white') +
    scale_fill_manual(na.value = 'white', breaks = 1:3, values = c("#a6cee3", "#1f78b4", "#b2df8a"), labels = c('Current', 'Future', 'Both')) +
    labs(fill = 'Most similar areas') +
    geom_polygon(data = shp, aes(x = long, y = lat, group = group), color = 'grey', fill = 'NA', position = ) +
    coord_equal() +
    #coord_equal(xlim = c(-89.35, -83.3), ylim =  c(lyr %>% extent() %>% as.vector() %>% .[3], 16.5)) +
    xlab('Longitude') + 
    ylab('Latitude') +
    theme_bw() + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'top') #north2(gg, x = .85, y = .85, symbol = 10) 
  ggplot2::ggsave(filename = outfile, plot = gg2, device = "jpeg", units = "in", width = 10, height = 6, dpi = 600)
  
}

# Summary table
map_summary <- function(current = T, crop = crop, stress = stress){
  
  if(current){
    lyr    <- prt_map
    period <- "Current"
  } else {
    lyr    <- ftr_map
    period <- "Future"
  }
  
  tmp <- lyr %>% raster::as.data.frame(xy = T)
  tmp <- tmp[complete.cases(tmp),]
  
  hrvs_clcl <- hrvs_area[][raster::cellFromXY(object = hrvs_area, xy = tmp %>% dplyr::select(x, y))] %>% sum(na.rm = T)
  ylds_clcl <- ylds_area[][raster::cellFromXY(object = ylds_area, xy = tmp %>% dplyr::select(x, y))] %>% sum(na.rm = T)
  prdc_clcl <- prdc_area[][raster::cellFromXY(object = prdc_area, xy = tmp %>% dplyr::select(x, y))] %>% sum(na.rm = T)
  
  data.frame(Crop = crop,
             Stress = stress,
             Period = period,
             Harvested_area = hrvs_clcl,
             Percent_harvested_area = paste0(round(hrvs_clcl/hrvs_area[] %>% sum(na.rm = T)*100, 2), "%"),
             Production = prdc_clcl,
             Percent_production = paste0(round(prdc_clcl/prdc_area[] %>% sum(na.rm = T)*100, 2), "%"),
             Yield = ylds_clcl,
             Percent_yield = paste0(round(ylds_clcl/ylds_area[] %>% sum(na.rm = T)*100, 2), "%")) %>% return()
}

# Calculate potential areas change
chg_areas_summary <- function(crop     = crop,
                              stress   = stress,
                              c_path   = c_path,
                              f_path   = f_path,
                              shp_path = shp_path,
                              outname  = outname,
                              c_areas  = c_areas)
  {
  
  # Present raster
  prt_map <<- tbl_to_rst(path = c_path, file = F, df = NULL, current = T)
  # Future raster
  ftr_map <<- min_ftr_dist(all_gcm = T, f_path = f_path, gcm_path = NULL)
  ftr_map <<- tbl_to_rst(path = NULL, file = T, df = ftr_map, current = F)
  
  # Change areas raster
  chg_map <<- raster::calc(x = stack(prt_map, ftr_map), fun = function(x) sum(x, na.rm = T))
  chg_map[which(chg_map[] == 0)] <- NA
  
  crop_areas <<- raster::brick(c_areas, lvar = 4)
  # ncdf4::nc_open("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_crop_presence/Bean/database/bean_monfreda.nc")
  hrvs_area <<- crop_areas[[5]] # Harvested area: Total hectares harvested per gridcell
  hrvs_area[which(hrvs_area[] == 0)] <- NA
  ylds_area <<- crop_areas[[2]] # Yield: Metric tons per hectare
  ylds_area[which(is.na(hrvs_area[]))] <- NA
  prdc_area <<- crop_areas[[6]] # Production: Metric tons
  prdc_area[which(is.na(hrvs_area[]))] <- NA
  
  # Plotting potential areas
  plot_maps(lyr = chg_map, shp_path = shp_path, outfile = outfile)
  
  # Summary table
  df <- lapply(c(T, F), function(x) map_summary(current = x, crop = crop, stress = stress)) %>% do.call(rbind, .)
  write.csv(df, file = "./summary_table_crop_stress.csv", row.names = F)
  return(df)
  
}

# ======================================================================= #
# Bean heat
# ======================================================================= #
crop     <- "Bean"
stress   <- "Heat"
c_path   <- paste0(root, "/CWR_pre-breeding/Results/Bean/Crop_indices/dist_levels_high_medium_low_ciat_palmira_corpoica.rds")
f_path   <- paste0(root, "/CWR_pre-breeding/Results/Bean/_future/rcp85")
shp_path <- paste0(root, "/Sustainable_Food_System/SFS_indicators/Input_data/world_shape/all_countries_edited.shp")
outfile  <- "./mapa_final.jpeg"
c_areas  <- paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Bean/database/bean_monfreda.nc")
chg_areas_summary(crop     = crop,
                  stress   = stress,
                  c_path   = c_path,
                  f_path   = f_path,
                  shp_path = shp_path,
                  outname  = outname,
                  c_areas  = c_areas)

# ======================================================================= #
# Sunflower heat
# ======================================================================= #
crop     <- "Sunflower"
stress   <- "Heat"
c_path   <- paste0(root, "/CWR_pre-breeding/Results/Sunflower/Crop_indices/dist_levels_high_medium_low_indian_narkhaoda.rds")
# f_path   <- paste0(root, "/CWR_pre-breeding/Results/Bean/_future/rcp85/gcm5/Crop_index")
f_path   <- paste0(root, "/CWR_pre-breeding/Results/Sunflower/_future/rcp85")
shp_path <- paste0(root, "/Sustainable_Food_System/SFS_indicators/Input_data/world_shape/all_countries_edited.shp")
outfile  <- "./sunflower_heat_map.jpeg"
chg_areas_summary(crop     = crop,
                  stress   = stress,
                  c_path   = c_path,
                  f_path   = f_path,
                  shp_path = shp_path,
                  outname  = outname,
                  c_areas  = c_areas)

# ======================================================================= #
# Carrot heat
# ======================================================================= #
crop     <- "Carrot"
stress   <- "Heat"
c_path   <- paste0(root, "/CWR_pre-breeding/Results/Carrot/Crop_indices/dist_levels_high_medium_low_thal_bhakar.rds")
f_path   <- paste0(root, "/CWR_pre-breeding/Results/Carrot/_future/heat/rcp85")
shp_path <- paste0(root, "/Sustainable_Food_System/SFS_indicators/Input_data/world_shape/all_countries_edited.shp")
outfile  <- "./carrot_heat_map.jpeg"
chg_areas_summary(crop     = crop,
                  stress   = stress,
                  c_path   = c_path,
                  f_path   = f_path,
                  shp_path = shp_path,
                  outname  = outname,
                  c_areas  = c_areas)

f_path   <- paste0(root, "/CWR_pre-breeding/Results/Carrot/_future/drought/rcp85")
