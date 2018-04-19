
suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})
suppressMessages(if(!require(FactoMineR)){install.packages('FactoMineR'); library(FactoMineR)} else {library(FactoMineR)})
suppressMessages(if(!require(FactoClass)){install.packages('FactoClass'); library(FactoClass)} else {library(FactoClass)})
suppressMessages(if(!require(ade4)){install.packages('ade4'); library(ade4)} else {library(ade4)})
suppressMessages(if(!require(xtable)){install.packages('xtable'); library(xtable)} else {library(xtable)})
suppressMessages(if(!require(ggdendro)){install.packages('ggdendro'); library(ggdendro)} else {library(ggdendro)})
suppressMessages(if(!require(compiler)){install.packages('compiler'); library(compiler)} else {library(compiler)})
suppressMessages(if(!require(ggthemes)){install.packages('ggthemes'); library(ggthemes)} else {library(ggthemes)})
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
suppressMessages(if(!require(googlesheets)){install.packages("googlesheets");library(googlesheets)}else{library(googlesheets)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})
suppressMessages(if(!require(caTools)){install.packages("caTools");library(caTools)}else{library(caTools)})
library(parallel)

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



index_carrot_drought_future <- function(continent = "Europa", rcp = "rcp85" ,gcm= "gcm1"){
  output <- paste0(root, "/CWR_pre-breeding/Results/Carrot/_future/drought/",rcp,"/",gcm,"/Crop_index/Carrot_index_drought_", tolower(continent), ".rds")
  if(!file.exists(output)){
    
    # Load climate data
    cat(">>> Starting process for carrot drought in", continent, "continent\n\n")
    cat(">>> Loading climate data ...\n")
    prec <- readRDS(paste0(root, '/CWR_pre-breeding/Input_data/_future_climate/',rcp,'/',gcm ,'/chirps/prec_filtered_', tolower(continent), '.rds'))
    
    prec$bean_coordinates <- NULL
    
    # Load crop cycle
    cat(">>> Loading crop cycle data ...\n")
    
    # Planting dates
    planting_rf_ggcmi <- raster::stack(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/carrot_planting.tif"))
    planting_rf_ggcmi <- planting_rf_ggcmi[[1]]
    
    # Harversting dates
    harvest_rf_ggcmi <- raster::brick(paste0(root, "/CWR_pre-breeding/Input_data/GGCMI-data/carrot_harvest.tif"))
    harvest_rf_ggcmi <- harvest_rf_ggcmi[[1]]
    
    
    # Extract important dates
    prec$Planting <- raster::extract(x = planting_rf_ggcmi, y = prec[,c("lon", "lat")])
    prec$Harvest  <- raster::extract(x = harvest_rf_ggcmi, y = prec[,c("lon", "lat")])
    prec$Duration <- ifelse(test = prec$Planting < prec$Harvest, yes = "One year", no = "Two years")
    
    
    # Restricting study area to crop area
    cat(">>> Restricting study area to crop area ...\n")
    crop_area  <- readRDS(paste0(root, "/CWR_pre-breeding/Input_data/_crop_presence/Carrot/database/area_base.rds"))
    prec <- dplyr::filter(prec, prec$cellID %in% crop_area$cellID)
    prec <- prec[!is.na(prec$cellID),]
    
    prec$Harvest[which(prec$Harvest == 0)] <- NA
    prec$Planting[which(prec$Planting == 0)] <- NA
    prec$Planting[which(prec$Harvest == -99)] <- NA
    prec<- na.omit(prec)
    prec$Harvest <- round(prec$Harvest)
    prec$Planting <- round(prec$Planting)
    
    
    test<- prec[which(prec$Duration =="Two years"), ]
    
    library(doSNOW)
    library(foreach)
    library(parallel)
    library(doParallel)
    
    cores<- detectCores()
    cl<- makeCluster(cores-18)
    registerDoParallel(cl) 
    
    system.time(indexes_drought <- foreach(i=1:nrow(prec)) %dopar% { ### 
      mylist <- list()
      cat(paste0("Processed pixel:", i, "\n"))
      # Parameters
      duration <- prec$Duration[i]
      start <- prec$Planting[i]
      end <- prec$Harvest[i]
      
      # Just one pixel
      time.serie <- prec[i, 1:(ncol(prec)-3)]
      
      if(duration== "One year"){
        suppressMessages(library(tidyverse))
        suppressMessages(library(compiler))
        
        ##NOTA :  X = Tmax , Y= Tmin , Z=Tmean
        
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% group_by(Year) %>% dplyr::filter(Yday >= start & Yday <= end)
        
        # TOTRAIN: Total precipitation
        
        totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
        totrain <- totrain %>% as.data.frame
        names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
        
        # prec_optimal 
        
        # 1.  350<prec<400
        
        prec_optimal <- totrain
        prec_optimal$con <- NA
        prec_optimal$con <- ifelse(prec_optimal$Value> 300  && prec_optimal$Value < 400, yes = 1, no = 0)
        prec_optimal <- data.frame(Year = prec_optimal$Year, Value= prec_optimal$con, Variable= "prec_optimal")
        
        # 2.   lack_prec <  250 
        # Drought
        
        lack_prec <- totrain
        lack_prec$con <- NA
        lack_prec$con <- ifelse(test = lack_prec$Value < 250 , yes = 1, no = 0)
        lack_prec <- data.frame(Year = lack_prec$Year, Value= lack_prec$con, Variable= "lack_prec")
        
        # 3. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
        dr_stress <- function(PREC, p_thresh = 1){
          runs <- rle(PREC < p_thresh)
          cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
          return(cons_days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
        cdd <- cdd %>% as.data.frame
        names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
        
        # 4. P_10: 95th percentile of daily precipitation (Erosion risk)
        p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .1, na.rm = TRUE))
        p_95 <- p_95 %>% as.data.frame
        names(p_95)[2] <- "Value"; p_95$Variable <- "P_10"
        
        # Putting all results together
        results <- data.frame(cellID = unique(X$cellID), rbind(prec_optimal, lack_prec, cdd, p_95))
        
      } else {
        
        suppressMessages(library(tidyverse))
        suppressMessages(library(compiler))
        
        X <- time.serie
        X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
        X$Year <- lubridate::year(as.Date(X$Date))
        X$Yday <- lubridate::yday(as.Date(X$Date))
        X <- X %>% filter(Yday %in% c(start:365, 1:end))
        X <- X[-(1:(end)),]
        X <- X[-((nrow(X)-(365-start)): nrow(X)),]
        X$condition <- NA
        for(j in 1:nrow(X)){
          X$condition[j] <- X$Yday[j+1] - X$Yday[j]
        }
        c <- c(unique(X$condition))
        X$condition[which(is.na(X$condition))] <- c[3]
        chngs <- which(X$condition == c[3])
        for(j in 1:length(chngs)){
          if(j == 1){
            X$condition2[1:chngs[j]] <- j
          } else {
            X$condition2[(chngs[j-1]+1):(chngs[j])] <- j
          }
        }
        X$condition2[is.na(X$condition2)] <- length(chngs)
        X$condition <- NULL
        names(X)[ncol(X)] <- "condition"
        rownames(X)<- 1:nrow(X)
        
        corte <- c(1:29)
        lista <- lapply(1:length(corte), function(i){
          tablas <- filter (X, condition== corte[i]) 
          return(tablas)  
        })
        
        X <- do.call(rbind,lista)
        X <- na.omit(X)
        
        # TOTRAIN: Total precipitation
        
        totrain <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
        totrain <- totrain %>% as.data.frame
        names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
        
        # prec_optimal 
        
        # 1.  350<prec<400
        
        prec_optimal <- totrain
        prec_optimal$con <- NA
        prec_optimal$con <- ifelse(prec_optimal$Value>300  && prec_optimal$Value < 400, yes = 1, no = 0)
        prec_optimal <- data.frame(Year = prec_optimal$Year, Value= prec_optimal$con, Variable= "prec_optimal")
        
        # 2.   lack_prec <  250 
        # Drought
        
        lack_prec <- totrain
        lack_prec$con <- NA
        lack_prec$con <- ifelse(test = lack_prec$Value < 250 , yes = 1, no = 0)
        lack_prec <- data.frame(Year = lack_prec$Year, Value= lack_prec$con, Variable= "lack_prec")
        
        # 3. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
        dr_stress <- function(PREC, p_thresh = 1){
          runs <- rle(PREC < p_thresh)
          cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
          return(cons_days)
        }
        dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
        cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
        cdd <- cdd %>% as.data.frame
        names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
        
        # 4. P_10: 95th percentile of daily precipitation (Erosion risk)
        p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .1, na.rm = TRUE))
        p_95 <- p_95 %>% as.data.frame
        names(p_95)[2] <- "Value"; p_95$Variable <- "P_10"
        
        # Putting all results together
        results <- data.frame(cellID = unique(X$cellID), rbind(prec_optimal, lack_prec, cdd, p_95))
        
      }
      mylist[[i]] <- results 
    } )
    
    stopCluster(cl)
    tabla <- do.call(rbind, indexes_drought)
    saveRDS(tabla, output)
    cat(">>> Results saved successfully ...\n")
    
    return(cat("Process done\n"))
    
  } else {
    cat(">>> Agroclimatic indices have been already calculated ...\n")
  }
  
}
system.time(index_carrot_drought_future(continent = "Africa", rcp = "rcp85" ,gcm= "gcm4"))