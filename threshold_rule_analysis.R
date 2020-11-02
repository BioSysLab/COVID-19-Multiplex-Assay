library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)
plate <- readRDS("data_rds/igg.rds")
type <- "igg"

# add labels

labels1 <- c(rep("Negative",78),rep("Positive",15))
labels2 <- c("Positive","Negative","Positive","Positive","Positive",
             "Negative","Negative","Positive","Positive","Negative",
             "Positive","Negative","Negative","Negative","Positive",
             "Positive","Negative","Negative","Positive","Negative",  #remember to set 18
             "Positive","Negative","Negative","Positive","Positive",
             "Negative","Negative","Positive","Positive","Negative",
             "Negative","Negative","Negative","Positive","Positive",
             "Positive","Positive","Negative","Negative","Negative",
             rep("Positive",53))
labels <- c(labels1,labels2)

# exclude the samples that are incorrect

plate$labels <- labels
excl <- c("48501","48524","48530")
plate <- plate[-which(plate$Sample%in%excl),]

# keep the columns that will be examined, depends on the concentration

columns <- c("REC31828.100_7.5UG","X40592.V08H_7.5UG","REC31812.100_30UG")
ncov <- plate[,c(1,ncol(plate),which(colnames(plate) %in% columns))]

# check if the colnames are in the correct order
###### CAREFUL HERE ##########

colnames(ncov) <- c("Sample","labels","S1","RBD","N")
ncov <- ncov[,c("Sample","labels","RBD","S1","N")]
conc <- c("0","0","7_ug","7_ug","7_ug")
conc <- c("0","0","15_ug","15_ug","15_ug")
conc <- c("0","0","30_ug","30_ug","30_ug")

###################################################

####### keep only the samples with known labels #############
ncov_all <- ncov
ncov <- ncov %>% filter(labels != "Blind samples")

cutoff <- calc_thresh(ncov,outlier = 1.5)


rbd_thresh <- seq(from = min(ncov$RBD),to = max(ncov$RBD),length.out = 2000)
s1_thresh <- seq(from = min(ncov$S1),to = max(ncov$S1), length.out = 2000)
n_thresh <- seq(from = min(ncov$N), to = max(ncov$N), length.out = 2000)

rbd_thresh <- seq(from = 2000,to = 5000,length.out = 1000)
s1_thresh <- seq(from = 7000,to = 8700, length.out = 1000)
n_thresh <- seq(from = 4800, to = 10600, length.out = 1000)


all_thresh <- list(0)
for (i in 1:length(rbd_thresh)) {
  df <- as.data.frame(matrix(0,nrow = 3, ncol = 2))
  colnames(df) <- c("predictor","roc_thresh")
  df$predictor <- c("RBD","S1","N")
  df$roc_thresh <- c(rbd_thresh[i],s1_thresh[i],n_thresh[i])
  all_thresh[[i]] <- df
}

results <- NULL
results <- foreach(thresh_n = all_thresh) %dopar% {
  eval_performance_roc_anal(ncov = ncov,roc_perf = thresh_n,
                       type = type,
                       conc = conc[3],
                       normal = "raw")
}

spec <- as.data.frame(matrix(0,nrow = length(results),ncol = 10))
sens <- as.data.frame(matrix(0,nrow = length(results),ncol = 10))

colnames(spec) <- c("rule","rbd","s1","n","rbd_s1","rbd_n","all","rbd_thresh","s1_thresh","n_thresh")
colnames(sens) <- c("rule","rbd","s1","n","rbd_s1","rbd_n","all","rbd_thresh","s1_thresh","n_thresh")

for (i in 1:length(results)) {
  df <- results[[i]]
  spec$rule[i] <- df$specificity[df$rule == "rule"]
  spec$rbd[i] <- df$specificity[df$rule == "RBD only"]
  spec$s1[i] <- df$specificity[df$rule == "S1 only"]
  spec$n[i] <- df$specificity[df$rule == "N only"]
  spec$rbd_s1[i] <- df$specificity[df$rule == "RBD + S1"]
  spec$rbd_n[i] <- df$specificity[df$rule == "RBD + N"]
  spec$all[i] <- df$specificity[df$rule == "all"]
  
  spec$rbd_thresh[i] <- df$thresh_rbd[df$rule == "RBD only"]
  spec$s1_thresh[i] <- df$thresh_s1[df$rule == "S1 only"]
  spec$n_thresh[i] <- df$thresh_n[df$rule == "N only"]
  
  sens$rule[i] <- df$sensitivity[df$rule == "rule"]
  sens$rbd[i] <- df$sensitivity[df$rule == "RBD only"]
  sens$s1[i] <- df$sensitivity[df$rule == "S1 only"]
  sens$n[i] <- df$sensitivity[df$rule == "N only"]
  sens$rbd_s1[i] <- df$sensitivity[df$rule == "RBD + S1"]
  sens$rbd_n[i] <- df$sensitivity[df$rule == "RBD + N"]
  sens$all[i] <- df$sensitivity[df$rule == "all"]
  
  sens$rbd_thresh[i] <- df$thresh_rbd[df$rule == "RBD only"]
  sens$s1_thresh[i] <- df$thresh_s1[df$rule == "S1 only"]
  sens$n_thresh[i] <- df$thresh_n[df$rule == "N only"]
  
}


plot(1:500,spec$rule,ylim = c(0.3,1.05),"o",col = "black",xlab = "id",
     ylab = "specificity" ,lwd = 1,lty = 5,pch = 1,,xaxs = "i",yaxs="i",cex = 0.6)

lines(1:500,spec$rbd,col = "#e41a1c",lwd = 1,type = "o",lty = 5,pch = 2,cex = 0.6)
lines(1:500,spec$s1,col = "#4daf4a",lwd = 1, type = "o",lty = 2 , pch = 3,cex = 0.6)
lines(1:500,spec$n,col = "#377eb8",lwd = 1, type = "o",lty = 5,pch = 4,cex = 0.6)
title("specificity",adj = 0)
legend("topleft", 
       legend = c("Rule", "N", "S1","RBD"),
       col = c('black', 
               '#377eb8',"#4daf4a","#e41a1c"),pch =c(1,4,3,2),
       pt.cex = 0.6, lwd = 0.8,
       cex = 0.6, 
       text.col = "black", 
       horiz = F )


plot(1:500,sens$rule,ylim = c(0.3,1.05),"o",col = "black",xlab = "id",
     ylab = "sensitivity" ,lwd = 1,lty = 5,pch = 1,,xaxs = "i",yaxs="i",cex = 0.6)

lines(1:500,sens$rbd,col = "#e41a1c",lwd = 1,type = "o",lty = 5,pch = 2,cex = 0.6)
lines(1:500,sens$s1,col = "#4daf4a",lwd = 1, type = "o",lty = 2 , pch = 3,cex = 0.6)
lines(1:500,sens$n,col = "#377eb8",lwd = 1, type = "o",lty = 5,pch = 4,cex = 0.6)
title("sensitivity",adj = 0)
legend("topleft", 
       legend = c("Rule", "N", "S1","RBD"),
       col = c('black', 
               '#377eb8',"#4daf4a","#e41a1c"),pch =c(1,4,3,2),
       pt.cex = 0.6, lwd = 0.8,
       cex = 0.6, 
       text.col = "black", 
       horiz = F )


library(MESS)

auc_rbd <- auc(x = (1-spec$rbd),y = sens$rbd,type = "linear")
auc_s1 <- auc(x = (1-spec$s1),y = sens$s1,type = "linear")
auc_n <- auc(x = (1-spec$n),y = sens$n,type = "linear")
auc_rule <- auc(x = (1-spec$rule),y = sens$rule,type = "linear")
auc_rbd_s1 <- auc(x = (1-spec$rbd_s1),y = sens$rbd_s1,type = "linear")
auc_rbd_n <- auc(x = (1-spec$rbd_n),y = sens$rbd_n,type = "linear")
auc_all <- auc(x = (1-spec$all),y = sens$all,type = "linear")

df_normal_auc <- c(auc_rbd,auc_s1,auc_n,auc_rule,auc_rbd_s1,auc_rbd_n,auc_all)
names(df_normal_auc) <- c("auc_rbd","auc_s1","auc_n","auc_rule","auc_rbd_s1","auc_rbd_n","auc_all")
write.csv(df_normal_auc,"roc_auc_all.csv")

w_spec <- 0.75
w_sens <- 0.25
plot(1:1000,w_spec*spec$rule+w_sens*sens$rule,ylim = c(0.9,1),"l",col = "black",xlab = "id",
     ylab = "specificity" ,lwd = 0.5,lty = 1,pch = 1,,xaxs = "i",yaxs="i",cex = 0.6)

lines(1:1000,w_spec*spec$rbd+w_sens*sens$rbd,col = "#e41a1c",lwd = 0.5,type = "l",lty = 5,pch = 2,cex = 0.6)
lines(1:1000,w_spec*spec$s1+w_sens*sens$s1,col = "#4daf4a",lwd = 0.5, type = "l",lty = 2 , pch = 3,cex = 0.6)
lines(1:1000,w_spec*spec$n+w_sens*sens$n,col = "#377eb8",lwd = 0.5, type = "l",lty = 5,pch = 4,cex = 0.6)
lines(1:1000,w_spec*spec$rbd_s1+w_sens*sens$rbd_s1,col = "purple",lwd = 0.5, type = "l",lty = 5,pch = 4,cex = 0.6)
lines(1:1000,w_spec*spec$rbd_n+w_sens*sens$rbd_n,col = "grey",lwd = 0.5, type = "l",lty = 5,pch = 4,cex = 0.6)
title(paste0("IGG ",w_spec,"*specificity+",w_sens,"*sensitivity"),adj = 0)
legend("bottomright", 
       legend = c("Rule", "N", "S1","RBD","RBD + S1","RBD + N"),
       col = c('black', 
               '#377eb8',"#4daf4a","#e41a1c","purple","grey"),pch =c(1,4,3,2),
       pt.cex = 0.6, lwd = 0.8,
       cex = 0.6, 
       text.col = "black", 
       horiz = F )
