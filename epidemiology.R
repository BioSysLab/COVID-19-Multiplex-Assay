library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
library(ggpubr)
library(RColorBrewer)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)
plate <- read.delim("data_exp2/total_epidem.txt")
type <- "total"

# label the samples

labels1 <- rep("Negative",nrow(plate)/2)
labels2 <- rep("Positive",nrow(plate)/2)
labels <- c(labels1,labels2)

# add to plate

plate$labels <- labels

# fix colnames of plate and ordering

colnames(plate) <- c("Sample","N","S1","RBD","labels")
ncov <- plate[,c("Sample","labels","RBD","S1","N")]

#set output dir
output_dir <- "paper_results/epidem"
conc <- "7ug"

perf_list <- eval_performance_neg_thresh(ncov = ncov,thresh = thresh,type = type,
                                         conc = conc,output_dir = output_dir,normal = "ratio",outlier_thresh = 0.3)

perf <- perf_list[[1]]
#write.csv(perf,paste0(output_dir,"/",type,"_performance.csv"))
predictions <- perf_list[[2]]

#save predictions

predictions <- left_join(ncov,predictions,by = "Sample") %>% dplyr::select(-labels.x) %>% dplyr::select(-labels.y)
predictions <- predictions %>% filter(RBD.y=="Positive"|s1=="Positive"|N.y=="Positive")

predictions <- predictions %>% select(Sample,RBD.x,S1,N.x,marginal_rule)

write.csv(predictions,paste0(output_dir,"/",type,"_predictions_4sd.csv"))

rbd_thresh <- seq(from = thresh$mean[1]+3*thresh$sd[1],to = thresh$mean[1]+5*thresh$sd[1],length.out = 200)
s1_thresh <- seq(from = thresh$mean[2]+3*thresh$sd[2],to = thresh$mean[2]+5*thresh$sd[2], length.out = 200)
n_thresh <- seq(from = thresh$mean[3]+3*thresh$sd[3],to = thresh$mean[3]+5*thresh$sd[3], length.out = 200)

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
  df <- results[[i]][[2]]
  #df <- df %>% select(rule,RBD,N,s1)
  id <- i
  prev_rule <- length(which(df$rule=='Positive'))/nrow(df)
  prev_RBD <- length(which(df$RBD=='Positive'))/nrow(df)
  prev_N <- length(which(df$N=='Positive'))/nrow(df)
  prev_S1 <- length(which(df$s1=='Positive'))/nrow(df)
  
  #prev_N_S1 <- length(which(df$n_and_s=='Positive'))/nrow(df)
  #prev_N_RBD <- length(which(df$n_rbd_rule=='Positive'))/nrow(df)
  #prev_RBD_S1 <- length(which(df$s1_rbd_rule=='Positive'))/nrow(df)
  
  prev_rule_N <- length(which(df$rule_n=='Positive'))/nrow(df)
  prev_rule_S <- length(which(df$rule_s=='Positive'))/nrow(df)
  prev_all <- length(which(df$all=='Positive'))/nrow(df)
  #prev_maj <- length(which(df$maj=='Positive'))/nrow(df)
  
  #prev_marg_pos <- length(which(df$marginal_positive=='Positive'))/nrow(df)
  #prev_marg_neg <- length(which(df$marginal_negative=='Positive'))/nrow(df)
  
  #df_spec <- df %>% dplyr::select(rule,specificity)
  #df_spec$id <- i
  #spec <- rbind(spec,df_spec)
  #df_sens <- df %>% dplyr::select(rule,sensitivity)
  #df_sens$id <- i
  #joined <- left_join(df_spec,df_sens,by = "rule")
  #,prev_rule_N,prev_rule_S,prev_all
  all_spec_sens <- rbind(all_spec_sens,c(id,prev_rule,prev_RBD,prev_N,prev_S1))
  #sens <- rbind(sens,df_sens)
}
all_spec_sens <- as.data.frame(all_spec_sens)
#,"N + RBD/S1","S1 + RBD/N","RBD + N + S1"
colnames(all_spec_sens) <- c("id","RBD & N|S1","RBD","N","S1")
#saveRDS(all_spec_sens,'prevalance_analysis_zoomed.rds')
#for epidemiology prevalance####
#breaks = c(0,100,200,300,800,1000),labels=c("0" = "Mean","200" = "Mean + 2SD","400" = "Mean + 4SD",
#"600" = "Mean + 6SD","800" = "Mean + 8SD","1000" = "Mean + 10SD")
#breaks = c(0,100,200,300),labels=c("0" = "Mean + 2SD","100" = "Mean + 3SD","200" = "Mean + 4SD",
 #                                  "300" = "Mean + 5SD")
#all_spec_sens <- all_spec_sens %>% select(prev_rule,prev_RBD,prev_N,prev_S1,id) %>% unique()
#colnames(all_spec_sens) <- c('RBD + N/S1','RBD','N','S1','id')
all_spec_sens <- gather(all_spec_sens,rule,prevalance,-id)
all_spec_sens$type <- "Combination"
all_spec_sens <- all_spec_sens %>% mutate(type = if_else((rule == "RBD")|(rule == "N")|(rule == "S1"),"Single",type))
titles <- "Total (IgG/IgA/IgM)"

# set factor order
#"S1 + RBD/N",,"N + RBD/S1","RBD + N + S1"
all_spec_sens$rule <- factor(all_spec_sens$rule,levels = c("RBD","RBD & N|S1","S1","N"))

cis <- NULL
cis <- lapply(all_spec_sens[,3]*1226, find_with_intervs_wilson,1255,1.96)

all_spec_sens$lower <- NULL
all_spec_sens$upper <- NULL
for (i in 1:nrow(all_spec_sens)) {
  all_spec_sens$lower[i] <- cis[[i]][1] 
  all_spec_sens$upper[i] <- cis[[i]][2]
}


library("ggsci")
png(filename = paste0(output_dir,"/","threshold_robustness/",type,"_prevalance_zoomed_ci.png"),width = 21,height = 14,units = "cm",res = 600)
#,300,400,500,600,700
#,"300" = "Mean + 6SD","400" = "Mean + 7SD","500" = "Mean + 8SD","600" = "Mean + 9SD","700" = "Mean + 10SD"
prevalance <-  ggplot(all_spec_sens, aes(x = id, y = prevalance*100)) +
  geom_line(aes(color = rule,linetype = type),position=position_dodge(width = 2),size = 0.8)+
  geom_ribbon(aes(ymin = lower*100,ymax = upper*100,fill = rule),alpha = 0.1)+
  scale_x_continuous(breaks = c(0,100,200),labels=c("0" = "Mean + 3SD","100" = "Mean + 4SD","200" = "Mean + 5SD")) +
  scale_y_continuous(breaks = seq(0,12,0.5),limits = c(0,12),expand = c(0,0))+
  xlab("Cut-off") + ylab("Estimated seroprevalence rate (%)")+theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=8))+
  theme(legend.text=element_text(size=8)) + 
  theme_minimal(base_size = 12)+
  #,"#33A02C","#E31A1C","#666666"
  labs(color = "SARS-CoV-2 antigen",fill = "95% confidence interval",linetype = "Decision type") + 
  scale_color_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#FB9A99"))+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#FB9A99"))
  
prevalance
dev.off()
##end of epidemiology prevalance


# start of plot for other corona

samples_n4 <- predictions %>% filter(N == "Positive") 
samples_n4 <- as.character(samples_n4$Sample)
samples_rule <- predictions %>% filter(rule == "Positive") 
samples_rule <- as.character(samples_rule$Sample)

samples_rbd <- predictions %>% filter(RBD == "Positive") 
samples_rbd <- as.character(samples_rbd$Sample)

samples_all <- predictions %>% filter(all == "Positive") 
samples_all <- as.character(samples_all$Sample)


samples_s <- predictions %>% filter(s1 == "Positive") 
samples_s <- as.character(samples_s$Sample)

raw <- read.csv("data_exp2/total_epidemio.csv")
raw <- raw[,c(1,5,6,7,8)]
colnames(raw) <- c("Sample","HKU S1","NL63 S1","229E S1","OC43 S1+S2")
raw_long <- reshape2::melt(raw)
library(ggrepel)

png(filename = paste0(output_dir,"/",type,"_N_pos_Rule_neg.png"),width = 9,height = 6,units = "in",res = 600)
ggplot(raw_long,aes(x=variable,y=value))+geom_violin(fill = "dodgerblue1")+
  geom_point(data=raw_long[which(raw_long$Sample %in% samples_n4),])+
  geom_label_repel(data =raw_long[which(raw_long$Sample %in% samples_n4),] ,aes(label=Sample),box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+ylab("Raw MFI") + xlab("") + ggtitle("N Positive and Rule Negative cases")
dev.off()

png(filename = paste0(output_dir,"/",type,"_N_neg_Rule_pos.png"),width = 9,height = 6,units = "in",res = 600)
ggplot(raw_long,aes(x=variable,y=value))+geom_violin(fill = "dodgerblue1")+
  geom_point(data=raw_long[which(raw_long$Sample %in% samples_rule),])+
  geom_label_repel(data =raw_long[which(raw_long$Sample %in% samples_rule),] ,aes(label=Sample),box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+ylab("Raw MFI") + xlab("") + ggtitle("Rule Positive and N Negative cases")
dev.off()

png(filename = paste0(output_dir,"/",type,"_S_pos_Rule_neg.png"),width = 9,height = 6,units = "in",res = 600)
ggplot(raw_long,aes(x=variable,y=value))+geom_violin(fill = "dodgerblue1")+
  geom_point(data=raw_long[which(raw_long$Sample %in% samples_s),])+
  geom_label_repel(data =raw_long[which(raw_long$Sample %in% samples_s),] ,aes(label=Sample),box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+ylab("Raw MFI") + xlab("") + ggtitle("Rule Positive and N Negative cases")
dev.off()

hku <- c(mean(raw$`HKU S1`[which(raw$Sample %in% samples_n4)]),
mean(raw$`HKU S1`[which(raw$Sample %in% samples_rule)]),
mean(raw$`HKU S1`[which(raw$Sample %in% samples_rbd)]),
mean(raw$`HKU S1`[which(raw$Sample %in% samples_s)]),
mean(raw$`HKU S1`))
hku <- as.data.frame(hku)
colnames(hku) <- "value"
hku$sd <- c(sd(raw$`HKU S1`[which(raw$Sample %in% samples_n4)]),
            sd(raw$`HKU S1`[which(raw$Sample %in% samples_rule)]),
            sd(raw$`HKU S1`[which(raw$Sample %in% samples_rbd)]),
            sd(raw$`HKU S1`[which(raw$Sample %in% samples_s)]),
            sd(raw$`HKU S1`))
hku$name <- "HKU S1"
hku$positive <- c("N","RBD + N/S1","RBD","S1","no subset")

nl63 <- c(mean(raw$`NL63 S1`[which(raw$Sample %in% samples_n4)]),
mean(raw$`NL63 S1`[which(raw$Sample %in% samples_rule)]),
mean(raw$`NL63 S1`[which(raw$Sample %in% samples_rbd)]),
mean(raw$`NL63 S1`[which(raw$Sample %in% samples_s)]),
mean(raw$`NL63 S1`))
nl63 <- as.data.frame(nl63)
colnames(nl63) <- "value"
nl63$sd <- c(sd(raw$`NL63 S1`[which(raw$Sample %in% samples_n4)]),
             sd(raw$`NL63 S1`[which(raw$Sample %in% samples_rule)]),
             sd(raw$`NL63 S1`[which(raw$Sample %in% samples_rbd)]),
             sd(raw$`NL63 S1`[which(raw$Sample %in% samples_s)]),
             sd(raw$`NL63 S1`))
nl63$name <- "NL63 S1"
nl63$positive <- c("N","RBD + N/S1","RBD","S1","no subset")

e <- c(mean(raw$`229E S1`[which(raw$Sample %in% samples_n4)]),
mean(raw$`229E S1`[which(raw$Sample %in% samples_rule)]),
mean(raw$`229E S1`[which(raw$Sample %in% samples_rbd)]),
mean(raw$`229E S1`[which(raw$Sample %in% samples_s)]),
mean(raw$`229E S1`))
e <- as.data.frame(e)
colnames(e) <- "value"
e$sd <- c(sd(raw$`229E S1`[which(raw$Sample %in% samples_n4)]),
          sd(raw$`229E S1`[which(raw$Sample %in% samples_rule)]),
          sd(raw$`229E S1`[which(raw$Sample %in% samples_rbd)]),
          sd(raw$`229E S1`[which(raw$Sample %in% samples_s)]),
          sd(raw$`229E S1`))
e$name <- "229E S1"
e$positive <- c("N","RBD + N/S1","RBD","S1","no subset")

oc43 <- c(mean(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_n4)]),
mean(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_rule)]),
mean(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_rbd)]),
mean(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_s)]),
mean(raw$`OC43 S1+S2`))
oc43 <- as.data.frame(oc43)
colnames(oc43) <- "value"
oc43$sd <- c(sd(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_n4)]),
             sd(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_rule)]),
             sd(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_rbd)]),
             sd(raw$`OC43 S1+S2`[which(raw$Sample %in% samples_s)]),
             sd(raw$`OC43 S1+S2`))

oc43$name <- "OC43 S1+S2"
oc43$positive <- c("N","RBD + N/S1","RBD","S1","no subset")

df <- bind_rows(e,oc43)
df$positive <- factor(df$positive,levels = c("no subset","RBD + N/S1","N","S1","RBD"))
png(filename = paste0(output_dir,"/",type,"meanmfi_othercov.png"),width = 9,height = 6,units = "in",res = 600)
ggplot(data=df, aes(x=name, y=value, fill=positive)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+xlab("")+ylab("Mean MFI")+labs(fill = "Antigen subset")+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9))
dev.off()
