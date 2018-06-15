# Validating waterlogging stress indices - CWR pre-breeding
# Author: H. Achicanoy
# CIAT, 2018

# Load libraries
library(tidyverse)
library(raster)

# Load drought stress indices by crop-stress and continent
# (It should vary by crop, stress and continent)
waterlogging_indices <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/Bean/Crop_indices/bean_waterlogging_crop_indices_europa.rds")
waterlogging_indices <- waterlogging_indices %>%
  dplyr::filter(Variable %in% c("Pr_ab_10", "Pr_ab_20", "Pr_ab_30", "Max_pr_day", "Totrain", "P5D"))

# Load crop cycle calendar by crop
# Planting date
planting_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "planting day")
planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
# Harvest day
harvest_rf_ggcmi <- raster::brick(paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/GGCMI-data/Pulses_rf_growing_season_dates_v1.25.nc4", sep = ""), varname = "harvest day")
harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]

# Identify pixel list where indices have been calculated
px_list <- waterlogging_indices$cellID %>% unique

# Load historical water balance for pixels where indices have been calculated
wb_files <- paste0("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_soils/Water_balance/Europa/", px_list, ".rds")
wb_files <- wb_files[which(file.exists(wb_files))]
px_list  <- px_list[which(file.exists(wb_files))]

# Sampling these files, choosing randomly 100 pixels for doing calculations
set.seed(1236)
wb_files_sub <- sample(wb_files, size = 100, replace = F)

# Calculating number of days within crop cycle with ERATIO < 0.5 (evidence of drought)
runoff <- wb_files_sub %>% purrr::map(.f = function(x){
  
  # Load water balance per pixel
  wb <- readRDS(x)
  wb$Year <- lubridate::year(as.Date(rownames(wb)))
  wb$Yday <- lubridate::yday(as.Date(rownames(wb)))
  
  # Extracting crop cycle
  start <- raster::extract(x = planting_rf_ggcmi, y = wb[,c("lon", "lat")] %>% unique)
  end   <- raster::extract(x = harvest_rf_ggcmi, y = wb[,c("lon", "lat")] %>% unique)
  
  # Extracting ERATIO variable
  runoff <- wb %>%
    group_by(Year) %>%
    dplyr::filter(Yday >= start & Yday <= end) %>%
    dplyr::arrange(Yday) %>%
    dplyr::summarise(Runoff = sd(RUNOFF, na.rm = T)) %>%
    as.data.frame
  runoff$cellID <- unique(as.numeric(wb$cellID))
  
  return(runoff)
  
})

# Identify selected pixels
slct_pixels <- gsub(pattern = "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Input_data/_soils/Water_balance/Europa/",
                    replacement = "",
                    wb_files_sub)
slct_pixels <- gsub(pattern = ".rds",
                    replacement = "",
                    slct_pixels)

# Wide format of drought indices
waterlogging_indices_wide <- waterlogging_indices %>% tidyr::spread(key = Variable, value = Value)
waterlogging_indices_wide$Max_pr_day <- as.numeric(waterlogging_indices_wide$Max_pr_day)
waterlogging_indices_wide$P5D <- as.numeric(waterlogging_indices_wide$P5D)
waterlogging_indices_wide$Pr_ab_10 <- as.numeric(waterlogging_indices_wide$Pr_ab_10)
waterlogging_indices_wide$Pr_ab_20 <- as.numeric(waterlogging_indices_wide$Pr_ab_20)
waterlogging_indices_wide$Pr_ab_30 <- as.numeric(waterlogging_indices_wide$Pr_ab_30)
waterlogging_indices_wide$Totrain <- as.numeric(waterlogging_indices_wide$Totrain)

# Calculating correlation between tentative drought indices vs ERATIO < 0.5 count
correlations <- lapply(1:length(slct_pixels), function(i){
  
  df   <- waterlogging_indices_wide %>% dplyr::filter(cellID == slct_pixels[i])
  df$cellID <- as.character(df$cellID)
  df$Year <- as.numeric(df$Year)
  runoff[[i]]$cellID <- as.character(runoff[[i]]$cellID)
  runoff[[i]]$Year <- as.numeric(runoff[[i]]$Year)
  m_df <- dplyr::left_join(x = df, y = runoff[[i]], by = c("cellID", "Year"))
  # summ <- data.frame(Variable = c("Pr_ab_10", "Pr_ab_20", "Pr_ab_30", "Max_pr_day", "Totrain", "P5D"),
  #                    Correlation = c(cor(m_df$Pr_ab_10,   m_df$Runoff, method = "pearson"),
  #                                    cor(m_df$Pr_ab_20,   m_df$Runoff, method = "pearson"),
  #                                    cor(m_df$Pr_ab_30,   m_df$Runoff, method = "pearson"),
  #                                    cor(m_df$Max_pr_day, m_df$Runoff, method = "pearson"),
  #                                    cor(m_df$Totrain,    m_df$Runoff, method = "pearson"),
  #                                    cor(m_df$P5D,        m_df$Runoff, method = "pearson")),
  #                    cellID = slct_pixels[i])
  
  return(m_df)
  
})
correlations <- do.call(rbind, correlations)
models_results <- correlations %>% dplyr::group_by(cellID) %>% do(broom::tidy(lm(Runoff ~ P5D + Totrain, data = .[,-c(1:2)])))
models_results %>% broom::tidy(fit)
models_results <- models_results[complete.cases(models_results),]
models_results <- models_results[models_results$p.value <= 0.05,]

table(models_results$estimate > 0, models_results$p.value <= 0.05, models_results$term)

models_results %>% dplyr::filter(term != "(Intercept)") %>% ggplot(aes(x = term, y = estimate)) + geom_boxplot()

correlations %>% ggplot(aes(y = Correlation, x = Variable)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, colour = "red") +
  geom_hline(yintercept = 0.3, colour = "blue")

correlations %>% modelr::fit_with()
