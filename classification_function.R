# Classification function
# Implemented by: H. Achicanoy
# CIAT, 2018

suppressMessages(library(pacman))
pacman::p_load(tidyverse, caret, randomForest, ISLR, nnet, caretEnsemble, ranger, ModelMetrics, rminer)

# sampling_mthd <- c("none", "down", "up")

## Run test
# df              <- all_data3
# standardize_all <- T
# sampling_mthd   <- "none"
# omit_correlated <- T
# top_variables   <- 5
# dscr_fun(df = df, standardize_all = standardize_all, sampling_mthd = sampling_mthd, omit_correlated = omit_correlated, top_variables = top_variables)

classification_fun <- function(df = all_data3,
                               standardize_all = T,
                               sampling_mthd = "none",
                               omit_correlated = T,
                               top_variables = 5,
                               ...)
{
  
  df_org <- df
  
  cat(">>> Pre-processing dataset ...\n")
  set.seed(1234)
  if(standardize_all){
    pp_rcp <- caret::preProcess(df_org[,-1], method = c("center", "scale", "nzv"))
    df <- predict(pp_rcp, df_org[,-1])
    df <- cbind(df_org[,1], df)
    names(df)[1] <- "Y"
  } else {
    pp_rcp <- preProcess(df_org[,-1], method = "nzv")
    df <- predict(pp_rcp, df_org[,-1])
    df <- cbind(df_org[,1], df)
    names(df)[1] <- "Y"
  }
  rm(df_org)
  
  if(omit_correlated){
    colinearity <- function(data){
      
      y              <- data[,1]
      numeric        <- data[,sapply(data, is.numeric)]
      descrCor       <- cor(numeric)
      highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = 0.8)
      
      if(length(highlyCorDescr) > 0){    
        numeric <- numeric[,-highlyCorDescr]
      }
      
      data <- data.frame(Y = y, numeric)
      return(data)
      
    }
    df <- colinearity(data = df)
  }
  
  cat(">>> Creating data partitions ...\n")
  inTrain    <- createDataPartition(y = df[,1], p = 0.7, list = FALSE)
  training   <- df[inTrain,]
  if(sampling_mthd == "down"){
    
    set.seed(1234)
    training <- caret::downSample(x = training[,-1],
                                  y = training[,1])
    training <- data.frame(Y = training$Class, training[,-which(colnames(training) == "Class")])
    
  } else {
    if(sampling_mthd == "up"){
      
      set.seed(1234)
      training <- caret::upSample(x = training[,-1],
                                  y = training[,1])
      training <- data.frame(Y = training$Class, training[,-which(colnames(training) == "Class")])
      
    } else {
      if(sampling_mthd == "none"){
        
        training <- training
        
      }
    }
  }
  testing    <- df[-inTrain,]
  
  cat(">>> Setting training parameters ...\n")
  control_prmt <- trainControl(method          = "LGOCV",
                               p               = 0.7,
                               number          = 10,
                               # index           = createResample(training$y, 10),
                               savePredictions = "final",
                               verboseIter     = T)
  
  cat(">>> Fitting 5 models: random forest, SVM, knn, FDA, avNNet ...\n")
  model_list <- caretEnsemble::caretList(
    Y ~ .,
    data       = training,
    trControl  = control_prmt,
    tuneList   = list(ranger = caretModelSpec(method = "ranger", importance = "impurity")),
    methodList = c("svmRadial", "knn", "bagFDA", "avNNet")
  )
  
  cat(">>> Calculating correlation and accuracy metrics over resamples ...\n")
  # It is expected high accuracy and un-correlation between them
  results <- list(model_correlations = modelCor(resamples(model_list)),
                  model_accuracies   = summary(resamples(model_list)))
  
  cat(">>> Calculating contingency table and metrics for training data ...\n")
  training_preds <- predict(model_list, newdata = training) %>% data.frame
  training_preds$ensemble <- apply(training_preds, 1, function(x){tt <- table(x); return(names(tt[which.max(tt)]))}) %>% factor
  results$Training_predictions <- training_preds
  
  results$Training_CM <- training_preds %>%
    purrr::map(function(x) suppressWarnings(caret::confusionMatrix(data = x, reference = training$Y)))
  
  results$Training_MCC <- training_preds %>%
    purrr::map(function(x) suppressWarnings(mltools::mcc(preds = x, actuals = training$Y))) %>%
    unlist
  
  cat(">>> Calculating contingency table and metrics for testing data ...\n")
  testing_preds <- predict(model_list, newdata = testing) %>% data.frame
  testing_preds$ensemble <- apply(testing_preds, 1, function(x){tt <- table(x); return(names(tt[which.max(tt)]))}) %>% factor
  results$Testing_predictions <- testing_preds
  
  results$Testing_CM <- testing_preds %>%
    purrr::map(function(x) suppressWarnings(caret::confusionMatrix(data = x, reference = testing$Y)))
  
  results$Testing_MCC <- testing_preds %>%
    purrr::map(function(x) suppressWarnings(mltools::mcc(preds = x, actuals = testing$Y))) %>%
    unlist
  
  cat(">>> Calculating variable importance for individual models ...\n")
  impVar_list <- lapply(1:length(model_list), function(i){
    vImportance <- caret::varImp(object = model_list[[i]])
    impVar <- data.frame(impVar = rownames(vImportance$importance)[1:top_variables])
    return(impVar)
  })
  impVar_list <- do.call(cbind, impVar_list)
  colnames(impVar_list) <- names(model_list)
  results$Important_variables <- impVar_list
  
  return(results)
  
}
