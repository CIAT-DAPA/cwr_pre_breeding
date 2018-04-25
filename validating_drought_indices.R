# Validating drought stress indices - CWR pre-breeding
# Author: H. Achicanoy
# CIAT, 2018

# Load libraries
library(tidyverse)
library(raster)

# Load drought stress indices by crop-stress and continent
# (It should vary by crop, stress and continent)
drought_indices <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/Lentil/Index_drought/lentil_index_drought_europa.rds")

# Load crop cycle calendar by crop
# Planting date
planting_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harvest day
harvest_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

# Identify pixel list where indices have been calculated
px_list <- drought_indices$cellID %>% unique

# Load historical water balance for pixels where indices have been calculated
wb_files <- paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_soils/Water_balance/Europa/", px_list, ".rds")
wb_files <- wb_files[which(file.exists(wb_files))]
px_list  <- px_list[which(file.exists(wb_files))]

# Sampling these files, choosing randomly 100 pixels for doing calculations
set.seed(1234)
wb_files_sub <- sample(wb_files, size = 100, replace = F)

# Calculating number of days within crop cycle with ERATIO < 0.5 (evidence of drought)
eratio <- wb_files_sub %>% purrr::map(.f = function(x){
  
  # Load water balance per pixel
  wb <- readRDS(x)
  wb$Year <- lubridate::year(as.Date(rownames(wb)))
  wb$Yday <- lubridate::yday(as.Date(rownames(wb)))
  
  # Extracting crop cycle
  start <- raster::extract(x = planting_rf_ggcmi, y = wb[,c("lon", "lat")] %>% unique)
  end   <- raster::extract(x = harvest_rf_ggcmi, y = wb[,c("lon", "lat")] %>% unique)
  
  # Extracting ERATIO variable
  eratio <- wb %>%
    group_by(Year) %>%
    dplyr::filter(Yday >= start & Yday <= end) %>%
    dplyr::arrange(Yday) %>%
    dplyr::summarise(Eratio = sum(ERATIO < 0.5, na.rm = T)) %>%
    as.data.frame
  eratio$cellID <- unique(as.numeric(wb$cellID))
  
  return(eratio)
  
})

# Identify selected pixels
slct_pixels <- gsub(pattern = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_soils/Water_balance/Europa/",
                    replacement = "",
                    wb_files_sub)
slct_pixels <- gsub(pattern = ".rds",
                    replacement = "",
                    slct_pixels)

# Wide format of drought indices
drought_indices_wide <- drought_indices %>% tidyr::spread(key = Variable, value = Value)

# Calculating correlation between tentative drought indices vs ERATIO < 0.5 count
correlations <- lapply(1:length(slct_pixels), function(i){
  
  df   <- drought_indices_wide %>% dplyr::filter(cellID == slct_pixels[i])
  m_df <- dplyr::left_join(x = df, y = eratio[[i]], by = c("cellID", "Year"))
  summ <- data.frame(Variable = c("prec_optimal", "lack_prec", "CDD", "P_10"),
                     Correlation = c(cor(m_df$prec_optimal, m_df$Eratio, method = "pearson"),
                                     cor(m_df$lack_prec,    m_df$Eratio, method = "pearson"),
                                     cor(m_df$CDD,          m_df$Eratio, method = "pearson"),
                                     cor(m_df$P_10,         m_df$Eratio, method = "pearson")),
                     cellID = slct_pixels[i])
  return(summ)
})
correlations <- do.call(rbind, correlations)
correlations %>% ggplot(aes(x = Correlation, fill = Variable)) + geom_density(alpha = .1) + geom_vline(xintercept = 0, colour = "red")
