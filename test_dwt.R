library(parallelDist)

# Multidimensional case without standarize
tmp.mat <- matrix(c(1:40), ncol = 10)
sample.matrix.list <- list(tmp.mat[1:2,], tmp.mat[3:4,])

parDist(x = sample.matrix.list, method = "dtw")

# Multidimensional case with standarize
sample.matrix.list.std <- sample.matrix.list
sample.matrix.list.std[[1]] <- t(apply(sample.matrix.list.std[[1]], MARGIN = 1, FUN = scale))
sample.matrix.list.std[[2]] <- t(apply(sample.matrix.list.std[[2]], MARGIN = 1, FUN = scale))

parDist(x = sample.matrix.list.std, method = "dtw")

# Unidimensional case without standarize
sample.matrix.list2 <- sample.matrix.list
sample.matrix.list2[[1]] <- cbind(as.matrix(t(sample.matrix.list2[[1]][1,])), as.matrix(t(sample.matrix.list2[[1]][2,])))
sample.matrix.list2[[2]] <- cbind(as.matrix(t(sample.matrix.list2[[2]][1,])), as.matrix(t(sample.matrix.list2[[2]][2,])))

parDist(x = sample.matrix.list2, method = "dtw")

# Our case
library(tidyverse)
df <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/Bean/General_indices/bean_general_indices_oceania.rds")
df <- df %>% filter(Variable %in% c("TMEAN", "GDD_1", "GDD_2", "ND_t35", "TOTRAIN", "CDD", "P5D", "P_95"))

df2 <- df %>% filter(cellID == "627688")
df2 <- df2 %>% spread(key = Year, value = Value)

df3 <- df %>% filter(cellID == "709769")
df3 <- df3 %>% spread(key = Year, value = Value)

colnames(df2)[3:ncol(df2)] <- paste0("Y", colnames(df2)[3:ncol(df2)])
colnames(df3)[3:ncol(df3)] <- paste0("Y", colnames(df3)[3:ncol(df3)])

sample.matrix.list3 <- list(data.matrix(df2[,3:ncol(df2)]), data.matrix(df3[,3:ncol(df3)]))
sample.matrix.list3[[1]] <- t(apply(X = sample.matrix.list3[[1]], MARGIN = 1, FUN = scale))
sample.matrix.list3[[2]] <- t(apply(X = sample.matrix.list3[[2]], MARGIN = 1, FUN = scale))

parDist(x = sample.matrix.list3, method = "dtw")

for(i in 1:nrow(df2)){
  l <- list(matrix(data = as.numeric(df2[i, 3:ncol(df2)]), nrow = 1, byrow = T),
            matrix(data = as.numeric(df3[i, 3:ncol(df3)]), nrow = 1, byrow = T))
  print(parDist(x = l, method = "dtw"))
}
sum(c(0.08918634, 0.008196721, 11, 1, 501, 1))
