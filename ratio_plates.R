library(tidyverse)

#total plate1

plate1 <- read.csv("performance_b/ratios/igm/plate1_igm.csv")
plate2 <- read.csv("performance_b/ratios/igm/plate2_igm.csv")

columns <- c("Sample","REC31828.100_15UG","REC31828.100_7.5UG","X40592.V08H_7.5UG","REC31812.100_30UG","X40592.V08H_15UG",
             "REC31828.100_30UG","X40592.V08H_30UG")

plate1 <- plate1[,which(colnames(plate1) %in% columns)]
plate2 <- plate2[,which(colnames(plate2) %in% columns)]

normalize_plate <- function(plate1){
  for (i in 2:ncol(plate1)) {
    if (!grepl(pattern = "REC31812",x = colnames(plate1)[i])) {
      plate1[,i] <- plate1[,i]/plate1[nrow(plate1)-1,i]
    } else {
      plate1[,i] <- plate1[,i]/plate1[nrow(plate1),i]
    }
    
  }
  nn <- nrow(plate1)
  plate1 <- plate1[-c(nn,nn-1),]
  return(plate1)
}

plate1 <- normalize_plate(plate1)
plate2 <- normalize_plate(plate2)

plate <- bind_rows(plate1,plate2)

saveRDS(plate,"data_rds/ratios/ratio_igm.rds")
