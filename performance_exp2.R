library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)
plate <- read.csv("data_exp2/igg_old.csv")
type <- "igg"

# calculate mean of negs

negs <- plate %>% filter(grepl(pattern = "Negative",x = Sample))
negs <- apply(negs[,2:4],MARGIN = 2,mean)

# remove the negatives and calculate ratios

plate <- plate %>% filter(!grepl(pattern = "Negative",x = Sample))
plate$REC31812 <- plate$REC31812/negs[1]
plate$REC31828 <- plate$REC31828/negs[2]
plate$X40592.V08H <- plate$X40592.V08H/negs[3]

# label the samples

labels1 <- rep("Negative",78)
labels2 <- rep("Positive",77)
labels <- c(labels1,labels2)

# add to plate

plate$labels <- labels

# fix colnames of plate and ordering

colnames(plate) <- c("Sample","N","S1","RBD","labels")
ncov <- plate[,c("Sample","labels","RBD","S1","N")]
# roc analysis
# set weights to sensitivity for Youden point calculation
w <- c(0.5,0.5,0.5)

#set output dir
output_dir <- "performance_exp2/igg"
conc <- "7ug"
# plot raw data
plot_scatter_bar_no_thresh(ncov = ncov,
                 output_dir = paste0(output_dir,"/plots"),
                 type = type,conc = conc)

roc_perf <- roc_anal(ncov = ncov,n = 1000,output_dir = paste0(output_dir,"/roc_analysis"),conc = conc,type = type,w=w)





performance_roc <- eval_performance_roc_anal(ncov = ncov,roc_perf = roc_perf,
                                             type = type,
                                             conc = "7ug",
                                             output_dir = paste0(output_dir,"/performance_eval"),
                                             normal = "ratio_neg", outlier_thresh = 0.4)

predictions <- performance_roc[[2]]

#write.csv(performance_roc,"performance_exp2/iga/performance_igg.csv")


# # try all combinations
# 
# rbd_thresh <- seq(from = min(ncov$RBD),to = max(ncov$RBD),length.out = 100)
# s1_thresh <- seq(from = min(ncov$S1),to = max(ncov$S1), length.out = 100)
# n_thresh <- seq(from = min(ncov$N), to = max(ncov$N), length.out = 100)
# 
# combinations <- expand.grid(rbd_thresh,s1_thresh,n_thresh)
# 
# colnames(combinations) <- c("rbd","s","n")
# 
# all_thresh <- list(0)
# for (i in 1:nrow(combinations)) {
#   df <- as.data.frame(matrix(0,nrow = 3, ncol = 2))
#   colnames(df) <- c("predictor","roc_thresh")
#   df$predictor <- c("RBD","S1","N")
#   df$roc_thresh <- c(combinations$rbd[i],combinations$s[i],combinations$n[i])
#   all_thresh[[i]] <- df
# }
# 
# results <- NULL
# results <- foreach(thresh_n = all_thresh) %dopar% {
#   eval_performance_roc_anal(ncov = ncov,roc_perf = thresh_n,
#                             type = type,
#                             conc = "7ug",
#                             normal = "neg_thresh")
# }


ncov_long <- reshape2::melt(ncov)
gscatters <- NULL
markers <- c("RBD","S1","N")
size <- c(10,6,5)
for (i in 1:3) {
  ncov_filt <- ncov_long %>% filter(variable == markers[i])
  dotplot <- ggplot(ncov_filt, aes(x=labels, y=value, fill = labels)) + 
    geom_dotplot(binaxis='y', 
                 stackdir='center', dotsize = size[i],method = "histodot", binwidth = 0.1)+
    xlab("")+ylab("ratio to negative") + ggtitle(paste0(type,"-",markers[i]))
  gscatters[[i]] <- dotplot
}

library(ggpubr)
png(filename = paste0(output_dir,"/",type,"_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
ggarrange(gscatters[[1]],gscatters[[2]],gscatters[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

# plot only outliers of rbd

out_rbd <- predictions[which(predictions$RBD_flag1==1),]
ncov_long_out <- ncov_long[which(ncov_long$Sample %in% out_rbd$Sample),]
gscatters_out <- NULL
markers <- c("RBD","S1","N")
size <- c(2,2,2)
for (i in 1:3) {
  ncov_filt <- ncov_long_out %>% filter(variable == markers[i])
  thresh <- roc_perf$roc_thresh[i]
  up <- 1.4 * thresh
  low <- 0.6 * thresh
  print(thresh)
  dotplot <- ggplot(ncov_filt, aes(x=labels, y=value, fill = labels)) + 
    geom_dotplot(binaxis='y', 
                 stackdir='center', dotsize = size[i],method = "histodot", binwidth = 0.1)+
    geom_hline(yintercept=thresh)+
    geom_hline(yintercept=up, linetype = "dashed")+
    geom_hline(yintercept=low, linetype = "dashed")+
    geom_text(aes(label=Sample),hjust=0, vjust=0)+
    ylim(0,15)+
    xlab("")+ylab("ratio to negative") + ggtitle(paste0(type,"-",markers[i]))
  gscatters_out[[i]] <- dotplot
}
png(filename = paste0(output_dir,"/",type,"_outlier_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
ggarrange(gscatters_out[[1]],gscatters_out[[2]],gscatters_out[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

plot_outliers_dist <- function(ncov,markers,roc_perf,marker){
  library(ggrepel)
  plotlist <- NULL
  ncov_long <- reshape2::melt(ncov)
  id_pos <- which(ncov$labels=="Positive")
  id_neg <- which(ncov$labels == "Negative")
  id_marker <- which(colnames(ncov)==marker)
  #mean_pos <- mean(ncov[id_pos,id_marker])
  #sd_pos <- sd(ncov[id_pos,id_marker])
  #mean_neg <- mean(ncov[id_neg,id_marker])
  #sd_neg <- sd(ncov[id_neg,id_marker])
  box_pos <- boxplot(ncov[id_pos,id_marker])
  lower_pos <- box_pos$stats[2]
  upper_pos <- box_pos$stats[4]
  iqr_pos <- 0.5*(upper_pos - lower_pos)
  thresh_pos <- lower_pos - iqr_pos
  box_neg <- boxplot(ncov[id_neg,id_marker])
  upper_neg <- box_neg$stats[4]
  lower_neg <- box_neg$stats[2]
  iqr_neg <- box_neg$stats[4]- box_neg$stats[2]
  thresh_neg <- upper_neg + iqr_neg
  samples_int_df <- ncov_long %>% filter(variable == marker) %>% filter(value <= thresh_pos) %>% filter(value >= thresh_neg)
  samples <- unique(as.character(samples_int_df$Sample))
  ncov_filt <- ncov_long[which(ncov_long$Sample %in% samples),]
  for (i in 1:length(markers)) {
    thresh <- roc_perf$roc_thresh[i]
    ncov_filt_plot <- ncov_filt %>% filter(variable == markers[i])
    id_pos <- which(ncov$labels=="Positive")
    id_neg <- which(ncov$labels == "Negative")
    id_marker <- which(colnames(ncov)==markers[i])
    #mean_pos <- mean(ncov[id_pos,id_marker])
    #sd_pos <- sd(ncov[id_pos,id_marker])
    #mean_neg <- mean(ncov[id_neg,id_marker])
    #sd_neg <- sd(ncov[id_neg,id_marker])
    box_pos <- boxplot(ncov[id_pos,id_marker])
    lower_pos <- box_pos$stats[2]
    upper_pos <- box_pos$stats[4]
    iqr_pos <- 0.5*(upper_pos - lower_pos)
    thresh_pos <- lower_pos - iqr_pos
    box_neg <- boxplot(ncov[id_neg,id_marker])
    upper_neg <- box_neg$stats[4]
    lower_neg <- box_neg$stats[2]
    iqr_neg <- box_neg$stats[4]- box_neg$stats[2]
    thresh_neg <- upper_neg + iqr_neg
    dotplot <- ggplot(ncov_filt_plot, aes(x=labels, y=value, fill = labels)) + 
      geom_dotplot(binaxis='y', 
                   stackdir='center', dotsize = 6,method = "histodot", binwidth = 0.1)+
      geom_hline(yintercept=thresh)+
      geom_hline(yintercept=thresh_pos, linetype = "dashed",colour = "blue")+
      geom_hline(yintercept=thresh_neg, linetype = "dashed",colour = "red")+
      geom_label_repel(aes(label=Sample),hjust=0, vjust=0,size = 2)+
      ylim(0,40)+
      xlab("")+ylab("ratio to negative") + ggtitle(paste0(type,"-",markers[i]))
    plotlist[[i]] <- dotplot
  }
  return(plotlist)
}
plotlist <- plot_outliers_dist(ncov,markers,marker = "RBD",roc_perf)
png(filename = paste0(output_dir,"/",type,"_distr_outlier_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
# rbd_thresh <- seq(from = 2,to = 3,length.out = 1000)
# s1_thresh <- seq(from = 0.75,to = 1.25, length.out = 1000)
# n_thresh <- seq(from = 0.6, to = 1.1, length.out = 1000)
# 
# 
# all_thresh <- list(0)
# for (i in 1:length(rbd_thresh)) {
#   df <- as.data.frame(matrix(0,nrow = 3, ncol = 2))
#   colnames(df) <- c("predictor","roc_thresh")
#   df$predictor <- c("RBD","S1","N")
#   df$roc_thresh <- c(rbd_thresh[i],s1_thresh[i],n_thresh[i])
#   all_thresh[[i]] <- df
# }
# 
# results <- NULL
# results <- foreach(thresh_n = all_thresh) %dopar% {
#   eval_performance_roc_anal(ncov = ncov,roc_perf = thresh_n,
#                             type = type,
#                             conc = conc[3],
#                             normal = "raw")
# }
# 
# spec <- as.data.frame(matrix(0,nrow = length(results),ncol = 10))
# sens <- as.data.frame(matrix(0,nrow = length(results),ncol = 10))
# 
# colnames(spec) <- c("rule","rbd","s1","n","rbd_s1","rbd_n","all","rbd_thresh","s1_thresh","n_thresh")
# colnames(sens) <- c("rule","rbd","s1","n","rbd_s1","rbd_n","all","rbd_thresh","s1_thresh","n_thresh")
# 
# for (i in 1:length(results)) {
#   df <- results[[i]]
#   spec$rule[i] <- df$specificity[df$rule == "rule"]
#   spec$rbd[i] <- df$specificity[df$rule == "RBD only"]
#   spec$s1[i] <- df$specificity[df$rule == "S1 only"]
#   spec$n[i] <- df$specificity[df$rule == "N only"]
#   spec$rbd_s1[i] <- df$specificity[df$rule == "RBD + S1"]
#   spec$rbd_n[i] <- df$specificity[df$rule == "RBD + N"]
#   spec$all[i] <- df$specificity[df$rule == "all"]
#   
#   spec$rbd_thresh[i] <- df$thresh_rbd[df$rule == "RBD only"]
#   spec$s1_thresh[i] <- df$thresh_s1[df$rule == "S1 only"]
#   spec$n_thresh[i] <- df$thresh_n[df$rule == "N only"]
#   
#   sens$rule[i] <- df$sensitivity[df$rule == "rule"]
#   sens$rbd[i] <- df$sensitivity[df$rule == "RBD only"]
#   sens$s1[i] <- df$sensitivity[df$rule == "S1 only"]
#   sens$n[i] <- df$sensitivity[df$rule == "N only"]
#   sens$rbd_s1[i] <- df$sensitivity[df$rule == "RBD + S1"]
#   sens$rbd_n[i] <- df$sensitivity[df$rule == "RBD + N"]
#   sens$all[i] <- df$sensitivity[df$rule == "all"]
#   
#   sens$rbd_thresh[i] <- df$thresh_rbd[df$rule == "RBD only"]
#   sens$s1_thresh[i] <- df$thresh_s1[df$rule == "S1 only"]
#   sens$n_thresh[i] <- df$thresh_n[df$rule == "N only"]
#   
# }
# 
# plot(1:nrow(spec),spec$rule,ylim = c(0.3,1.05),"o",col = "black",xlab = "id",
#      ylab = "specificity" ,lwd = 1,lty = 5,pch = 1,,xaxs = "i",yaxs="i",cex = 0.6)
# 
# lines(1:nrow(spec),spec$rbd,col = "#e41a1c",lwd = 1,type = "o",lty = 5,pch = 2,cex = 0.6)
# lines(1:nrow(spec),spec$s1,col = "#4daf4a",lwd = 1, type = "o",lty = 2 , pch = 3,cex = 0.6)
# lines(1:nrow(spec),spec$n,col = "#377eb8",lwd = 1, type = "o",lty = 5,pch = 4,cex = 0.6)
# title("specificity",adj = 0)
# legend("bottomright", 
#        legend = c("Rule", "N", "S1","RBD"),
#        col = c('black', 
#                '#377eb8',"#4daf4a","#e41a1c"),pch =c(1,4,3,2),
#        pt.cex = 0.6, lwd = 0.8,
#        cex = 0.6, 
#        text.col = "black", 
#        horiz = F )
# 
# plot(1:nrow(sens),sens$rule,ylim = c(0.3,1.05),"o",col = "black",xlab = "id",
#      ylab = "sensitivity" ,lwd = 1,lty = 5,pch = 1,,xaxs = "i",yaxs="i",cex = 0.6)
# 
# lines(1:nrow(sens),sens$rbd,col = "#e41a1c",lwd = 1,type = "o",lty = 5,pch = 2,cex = 0.6)
# lines(1:nrow(sens),sens$s1,col = "#4daf4a",lwd = 1, type = "o",lty = 2 , pch = 3,cex = 0.6)
# lines(1:nrow(sens),sens$n,col = "#377eb8",lwd = 1, type = "o",lty = 5,pch = 4,cex = 0.6)
# title("sensitivity",adj = 0)
# legend("topleft", 
#        legend = c("Rule", "N", "S1","RBD"),
#        col = c('black', 
#                '#377eb8',"#4daf4a","#e41a1c"),pch =c(1,4,3,2),
#        pt.cex = 0.6, lwd = 0.8,
#        cex = 0.6, 
#        text.col = "black", 
#        horiz = F )
