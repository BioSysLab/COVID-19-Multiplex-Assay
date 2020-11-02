library(tidyverse)

plate1 <- read.csv("performance_b/plate1_iga.csv")
#plate1 <- plate1[,-1]
plate1$Sample <- as.character(plate1$Sample)

plate2 <- read.csv("performance_b/plate2_iga.csv")
plate2$Sample <- as.character(plate2$Sample)

labels1 <- c(rep("Negative",78),rep("Positive",15))
labels2 <- c(rep("Blind samples",40),rep("Positive",53))
labels <- c(labels1,labels2)

igg <- bind_rows(plate1,plate2)
igg$labels <- labels
#"REC31806.500_30UG"
ncov <- c("REC31828.100_30UG","X40592.V08H_30UG","REC31812.100_30UG")
ncov <- igg[,c(1,21,which(colnames(igg) %in% ncov))]

ncov_test <- ncov %>% filter(labels == "Blind samples")
colnames(ncov_test) <- c("Sample","labels","N","S1_old","RBD")

ncov_val <- ncov %>% filter(labels != "Blind samples")
#samples to exclude
excl <- c("48501","48524","48530")
ncov_val <- ncov_val[-which(ncov_val$Sample%in%excl),]
colnames(ncov_val) <- c("Sample","labels","N","S1_old","RBD")


# the parameter to consider outlier is important here dont forget#####
######################################################################
######################################################################
calc_thresh <- function(df){
  #first is always Samples and second is always labels
  n <- ncol(df)
  npred <- n-2
  limits <- matrix(0,nrow = npred,ncol = 4)
  limits <- as.data.frame(limits)
  colnames(limits) <- c("predictor","lower","thresh","upper")
  for (i in 1:npred) {
    meas <- as.numeric(df[,i+2])
    meas_neg <- meas[df$labels=="Negative"]
    m <- mean(meas_neg)
    sd <- sd(meas_neg)
    thresh <- m + 3*sd
    id <- which(meas[df$labels=="Negative"]>(1.5)*thresh)
    if (length(id) != 0) {
      neg <- which(df$labels=="Negative")
      neg <- neg[-id]
      thresh <- mean(meas[neg]) + 3* sd(meas[neg])
      m <- mean(meas[neg])
      sd <- sd(meas[neg])
    }
    limits[i,2] <- m + 2*sd
    limits[i,3] <- m + 3*sd
    limits[i,4] <- m + 4*sd
  }
  limits$predictor <- as.character(colnames(df)[3:n])
  return(limits)
}

limit <- calc_thresh(ncov_val)
npred <- ncol(ncov_val)-2

predictions <- NULL
predictions <- ncov_test[,c(1,2)]
predictions <- as.data.frame(predictions)
for (i in 1:npred) {
  meas <- as.numeric(ncov_test[,i+2])
  pred <- meas >= limit$thresh[i]
  pred <- pred + 0
  flag1 <- (meas>=limit$lower[i] & meas <= limit$upper[i])
  flag1 <- flag1 + 0
  predictions <- cbind(predictions,pred,flag1)
}

s <- seq(1,npred*2,2)
#naming the predictions
for (i in 1:length(s)) {
  
  colnames(predictions)[s[i]+2] <- paste0(colnames(ncov_test)[i+2],"_pred")
  colnames(predictions)[s[i]+3] <- paste0(colnames(ncov_test)[i+2],"_flag1")
  
}

predictions_rule <- predictions %>% mutate(pred_agr = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_old_pred == 1 )),
                                                          true = "Positive",false = "Negative"))

predictions_rule$flag2 <- 0
predictions_rule$flag2[predictions_rule$pred_agr == "Positive" & predictions_rule$RBD_flag1 ==1] <- 1


predictions_rbd <- predictions %>% mutate(pred_agr = if_else(RBD_pred==1,true = "Positive",false = "Negative"))

predictions_rule$pred_RBD <- predictions_rbd$pred_agr

which(predictions_rule$pred_agr != predictions_rule$pred_RBD)

write.csv(predictions_rule,"iga_results/test_predictions.csv")
