library(tidyverse)

plate1 <- read.csv("performance_b/plate1_igg_igm_iga.csv")
#plate1 <- plate1[,-1]
plate1$Sample <- as.character(plate1$Sample)

plate2 <- read.csv("performance_b/plate2_igg_igm_iga.csv")
plate2$Sample <- as.character(plate2$Sample)

labels1 <- c(rep("Negative",78),rep("Positive",15))
labels2 <- c(rep("Blind samples",40),rep("Positive",53))
labels <- c(labels1,labels2)
#labels2 <- c(rep("negative_gennimatas",66),rep("negative_patras",12),rep("positive_patras",15),rep("test",40),rep("positive_alexandras",53))


igg <- bind_rows(plate1,plate2)
igg$labels <- labels
#"REC31806.500_30UG"
ncov <- c("REC31828.100_30UG","X40592.V08H_30UG","REC31812.100_30UG")
ncov <- igg[,c(1,21,which(colnames(igg) %in% ncov))]
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
    id <- which(meas[df$labels=="Negative"]>(2.5)*thresh)
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
#results <- matrix(0,nrow = nrow(ncov_val), ncol = (2+2*npred))
#results <- as.data.frame(results)
#results$V1 <- as.character(ncov_val$Sample)
#results$V2 <- as.character(ncov_val$labels)
# make predictions from each predictor and flag the iffy ones
predictions <- NULL
predictions <- ncov_val[,c(1,2)]
predictions <- as.data.frame(predictions)
for (i in 1:npred) {
  meas <- as.numeric(ncov_val[,i+2])
  pred <- meas >= limit$thresh[i]
  pred <- pred + 0
  flag1 <- (meas>=limit$lower[i] & meas <= limit$upper[i])
  flag1 <- flag1 + 0
  predictions <- cbind(predictions,pred,flag1)
}

s <- seq(1,npred*2,2)
#naming the predictions
for (i in 1:length(s)) {
  
  colnames(predictions)[s[i]+2] <- paste0(colnames(ncov_val)[i+2],"_pred")
  colnames(predictions)[s[i]+3] <- paste0(colnames(ncov_val)[i+2],"_flag1")
  
}

# aggregate predictions

predictions1 <- predictions %>% mutate(pred_agr = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_old_pred == 1 )),
                                                         true = "Positive",false = "Negative"))
#df <- predictions1
eval_preds <- function(df){
  tp <- df %>% filter(pred_agr == "Positive") %>% filter(labels=="Positive")
  tp <- nrow(tp)
  
  fp <- df %>% filter(pred_agr == "Positive") %>% filter(labels=="Negative")
  fp <- nrow(fp)
  
  tn <- df %>% filter(pred_agr == "Negative") %>% filter(labels=="Negative")
  tn <- nrow(tn)
  
  fn <- df %>% filter(pred_agr == "Negative") %>% filter(labels=="Positive")
  fn <- nrow(fn)
  
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  ppv <- tp/(tp+fp)
  npv <- tn/(tn+fn)
  acc <- (tn+tp)/(tn+fn+tp+fp)
  # add confidence intervals
  return(c(sens, spec, ppv,npv, acc,fn,fp))
}

results1 <- eval_preds(predictions1)
predictions2 <- predictions %>% mutate(pred_agr = if_else(RBD_pred==1,true = "Positive",false = "Negative"))
results2 <- eval_preds(predictions2)

names(results1) <- c("sensitivity","specificity","ppv","npv","accuracy","fn","fp")
names(results2) <- c("sensitivity","specificity","ppv","npv","accuracy","fn","fp")

write.csv(results1,"iga_results/iga_metrics_rule_based.csv")
write.csv(results2,"iga_results/ig_all_metrics_RBD_only.csv")

predictions1$flag2 <- 0
predictions1$flag2[predictions1$pred_agr == "Positive" & predictions1$RBD_flag1 ==1] <- 1

write.csv(predictions1,"iga_results/iga_val_predictions.csv")

#predictions2 <- predictions
#predictions2$N_conf <- -1
#predictions2$N_conf[(predictions2$N_pred==1 & predictions2$N_flag1==1)] <- 0.5
#predictions2$N_conf[(predictions2$N_pred==0 & predictions2$N_flag1==1)] <- -0.5
#predictions2$N_conf[(predictions2$N_pred==1 & predictions2$N_flag1==0)] <- 1

#predictions2$S1_conf <- -1
#predictions2$S1_conf[(predictions2$S1_old_pred==1 & predictions2$S1_old_flag1==1)] <- 0.5
#predictions2$S1_conf[(predictions2$S1_old_pred==0 & predictions2$S1_old_flag1==1)] <- -0.5
#predictions2$S1_conf[(predictions2$S1_old_pred==1 & predictions2$S1_old_flag1==0)] <- 1


#predictions2$RBD_conf <- -1
#predictions2$RBD_conf[(predictions2$RBD_pred==1 & predictions2$RBD_flag1==1)] <- 0.5
#predictions2$RBD_conf[(predictions2$RBD_pred==0 & predictions2$RBD_flag1==1)] <- -0.5
#predictions2$RBD_conf[(predictions2$RBD_pred==1 & predictions2$RBD_flag1==0)] <- 1

#predictions2$pred_agr <- apply(predictions2[,9:11],1,sum)

#predictions2$pred_agr[predictions2$pred_agr>0] <- 1
#predictions2$pred_agr[predictions2$pred_agr<0] <- -1
#predictions2$pred_agr[predictions2$pred_agr==0] <- predictions2$RBD_pred[predictions2$pred_agr==0]

#predictions2 <- predictions2 %>% mutate(pred_agr = if_else(pred_agr == 1,
 #                                                         true = "Positive",false = "Negative"))


#results2 <- eval_preds(predictions2)

#first is always Samples and second is always labels

predictions_rule <- predictions
predictions_rule$pred_agr <- "empty"
for (i in 1:nrow(predictions_rule)) {
  if (predictions_rule$RBD_flag1[i]==0) {
    if (predictions_rule$RBD_pred[i]==1){
      predictions_rule$pred_agr[i] <- "Positive"
    }
    if (predictions_rule$RBD_pred[i]==0){
      predictions_rule$pred_agr[i] <- "Negative"
    }
  }
  if (predictions_rule$RBD_flag1[i]==1) {
    if (predictions_rule$RBD_pred[i]==1 & (predictions_rule$N_pred[i] == 1 | predictions_rule$S1_old_pred[i] == 1)){
      predictions_rule$pred_agr[i] <- "Positive"
    } else {
      predictions_rule$pred_agr[i] <- "Negative"
    }
  }
}
results3 <- eval_preds(predictions_rule)
