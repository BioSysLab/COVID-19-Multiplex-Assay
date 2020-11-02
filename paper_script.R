library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
library(ggpubr)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)
plate <- read.csv("data_exp2/total_old.csv")
type <- "total"

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

#set output dir
output_dir <- "paper_results/total"
conc <- "7ug"
threshold_sds <- "3_wilson"

thresh <- calc_neg_thresh_wo_outliers(ncov,2,thresh_sd=3)
#thresh <- calc_truncated_neg_thresh_wo_outliers(df = ncov,iqr = 3,thresh_sd = 3)

# plot raw data
plot_scatter_bar_neg_thresh(ncov = ncov,thresh,
                           output_dir = paste0(output_dir,"/plots"),
                           type = type,conc = conc)


# check the plotted thresh

ncov_long <- reshape2::melt(ncov)
plotlist <- NULL
markers <- c("RBD","S1","N")
size <- c(10,10,10)
for (i in 1:3) {
  ncov_filt <- ncov_long %>% filter(variable == markers[i])
  thresh_neg <-  thresh$thresh[i]
  dotplot <- ggplot(ncov_filt, aes(x=labels, y=value, fill = labels)) + 
    geom_dotplot(binaxis='y', 
                 stackdir='center', dotsize = size[i],method = "histodot", binwidth = 0.1)+
    geom_hline(yintercept=thresh_neg, linetype = "dashed",colour = "red")+
    ylim(0,max(ncov_long$value))+
    xlab("")+ylab("ratio to negative") + ggtitle(paste0(type,"-",markers[i]))
  plotlist[[i]] <- dotplot
}

png(filename = paste0(output_dir,"/",type,"_dotplots_w_neg_threshold.png"),width = 18,height = 6,units = "in",res = 600)
ggarrange(plotlist[[1]],plotlist[[2]],plotlist[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

perf_list <- eval_performance_neg_thresh(ncov = ncov,thresh = thresh,type = type,
                                    conc = conc,output_dir = output_dir,normal = "ratio",outlier_thresh = 0.3)

perf <- perf_list[[1]]
write.csv(perf,paste0(output_dir,"/",type,threshold_sds,"_performance.csv"))
predictions <- perf_list[[2]]


rbd_thresh <- seq(from = thresh$mean[1]+thresh$sd[1],to = thresh$mean[1]+5*thresh$sd[1],length.out = 400)
s1_thresh <- seq(from = thresh$mean[2]+thresh$sd[2],to = thresh$mean[2]+5*thresh$sd[2], length.out = 400)
n_thresh <- seq(from = thresh$mean[3]+thresh$sd[3],to = thresh$mean[3]+5*thresh$sd[3], length.out = 400)

# rbd_thresh <- seq(from = min(ncov$RBD),to = max(ncov$RBD),length.out = 2000)
# s1_thresh <- seq(from = min(ncov$S1),to = max(ncov$S1), length.out = 2000)
# n_thresh <- seq(from = min(ncov$N), to = max(ncov$N), length.out = 2000)

all_thresh <- list(0)
for (i in 1:length(rbd_thresh)) {
  df <- as.data.frame(matrix(0,nrow = 3, ncol = 2))
  colnames(df) <- c("predictor","thresh")
  df$predictor <- c("RBD","S1","N")
  df$thresh <- c(rbd_thresh[i],s1_thresh[i],n_thresh[i])
  all_thresh[[i]] <- df
}

results <- NULL
results <- foreach(thresh_n = all_thresh) %dopar% {
  eval_performance_neg_thresh(ncov = ncov,thresh = thresh_n,
                            type = type,
                            conc = conc,
                            normal = "ratio", outlier_thresh = 0.4)
}

all_spec_sens <- NULL
for (i in 1:length(results)) {
  df <- results[[i]][[1]]
  df_spec <- df %>% dplyr::select(rule,specificity,accuracy)
  df_spec$id <- i
  #spec <- rbind(spec,df_spec)
  df_sens <- df %>% dplyr::select(rule,sensitivity)
  #df_sens$id <- i
  joined <- left_join(df_spec,df_sens,by = "rule")
  all_spec_sens <- rbind(all_spec_sens,joined)
  #sens <- rbind(sens,df_sens)
}
#"RBD + S1","RBD + N","N + S","rule N","rule S","Majority",,"RBD + S1","all", "rule S"
rules <- c("rule","RBD only","S1 only","N only","rule S","rule N","all")
all_filtered <- all_spec_sens[which(all_spec_sens$rule %in% rules),]
all_filtered$sum <- (all_filtered$sensitivity+all_filtered$specificity)/2
all_filtered$one_spec <- 1 - all_filtered$specificity
all_filtered$rule[all_filtered$rule == "all"] <- "RBD & N & S1"
all_filtered$rule[all_filtered$rule == "rule"] <- "RBD & N|S1"
all_filtered$rule[all_filtered$rule == "RBD only"] <- "RBD"
all_filtered$rule[all_filtered$rule == "S1 only"] <- "S1"
all_filtered$rule[all_filtered$rule == "N only"] <- "N"
all_filtered$rule[all_filtered$rule == "rule S"] <- "S1 & RBD|N"
all_filtered$rule[all_filtered$rule == "rule N"] <- "N & RBD|S1"

all_filtered <- all_filtered %>% mutate(type = if_else(rule == "N" | rule == "S1" | rule == "RBD",true = "Single",false = "Combination"))

titles <- "Total (IgG/IgA/IgM)"
titles <- "IgG"
titles <- "IgA"
titles <- "IgM"

# set factor order
all_filtered$rule <- factor(all_filtered$rule,levels = c("RBD","RBD & N|S1","S1","S1 & RBD|N","N","N & RBD|S1","RBD & N & S1"))


png(filename = paste0(output_dir,"/","threshold_robustness/",type,"_specificity.png"),width = 9,height = 6,units = "in",res = 600)

#scale_x_continuous(breaks = c(0,100,200,300,400),labels=c("0" = "Mean + SD", "100" = "Mean + 2SD","200" = "Mean + 3SD",
#                                                          "300" = "Mean + 4SD","400" = "Mean + 5SD")) 

#scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200),labels=c("0" = "Mean + SD", "100" = "","200" = "",
#"300" = "Mean + 4SD","400" = "","500" = "","600" = "Mean + 7SD","700"="",
#"800" = "","900"="Mean + 10SD","1000" = "","1100"="","1200" = "Mean + 13SD"))

specificity <- ggplot(all_filtered, aes(x = id, y = specificity)) + 
    geom_line(aes(color = rule,linetype = type),position=position_dodge(width = 6),size = 0.8)+ylim(0.7,1)+
  scale_x_continuous(breaks = c(0,200,400),labels=c("0" = "Mean + SD", "200" = "Mean + 3SD",
                                                    "400" = "Mean + 5SD"))  + xlab("Cut-off") + ylab("Specificity")+
    theme(legend.text=element_text(size=10),axis.text.x = element_text(size = 10)) + labs(color = "SARS-CoV-2 antigen",linetype = "Decision type")+
  theme_minimal(base_size = 18) + scale_color_manual(values = c("#A6CEE3","#1F78B4",
                                                                "#B2DF8A","#33A02C",
                                                                "#FB9A99","#E31A1C",
                                                                "#666666"))
specificity
dev.off()

png(filename = paste0(output_dir,"/","threshold_robustness/",type,"_sensitivity.png"),width = 9,height = 6,units = "in",res = 600)
sensitivity <- ggplot(all_filtered, aes(x = id, y = sensitivity)) + 
    geom_line(aes(color = rule,linetype = type),position=position_dodge(width = 6),size = 0.8)+ylim(0.7,1)+
  scale_x_continuous(breaks = c(0,200,400),labels=c("0" = "Mean + SD", "200" = "Mean + 3SD",
                                                    "400" = "Mean + 5SD"))  + xlab("Cut-off") + xlab("Cut-off") + ylab("Sensitivity")+
  theme(legend.text=element_text(size=10),axis.text.x = element_text(size = 10)) + labs(color = "SARS-CoV-2 antigen",linetype = "Decision type")+
  theme_minimal(base_size = 18)+ scale_color_manual(values = c("#A6CEE3","#1F78B4",
                                                               "#B2DF8A","#33A02C",
                                                               "#FB9A99","#E31A1C",
                                                               "#666666"))
sensitivity
dev.off()

png(filename = paste0(output_dir,"/","threshold_robustness/",type,"_accuracy.png"),width = 9,height = 6,units = "in",res = 600)
acc_plot <- ggplot(all_filtered, aes(x = id, y = accuracy)) + 
  geom_line(aes(color = rule,linetype = type),position=position_dodge(width = 6),size = 0.8)+ylim(0.7,1)+
  scale_x_continuous(breaks = c(0,200,400),labels=c("0" = "Mean + SD", "200" = "Mean + 3SD",
                                                    "400" = "Mean + 5SD"))  + xlab("Cut-off") + xlab("Cut-off") + ylab("Specificity")+
  ylab("Accuracy")+
  theme(legend.text=element_text(size=10),axis.text.x = element_text(size = 10)) + labs(color = "SARS-CoV-2 antigen",linetype = "Decision type")+
  theme_minimal(base_size = 18)+ scale_color_manual(values = c("#A6CEE3","#1F78B4",
                                                               "#B2DF8A","#33A02C",
                                                               "#FB9A99","#E31A1C",
                                                               "#666666"))
acc_plot
dev.off()

png(filename = paste0(output_dir,"/","threshold_robustness/",type,"_combined_leg.png"),width = 40,height = 12,units = "cm",res = 600)
combined <- ggarrange(specificity,sensitivity,acc_plot, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
annotate_figure(combined)
dev.off()
