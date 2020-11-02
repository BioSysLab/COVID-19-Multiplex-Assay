library(tidyverse)
library(viridis)
library(reshape2)

# load plate

plate <- readRDS("data_rds/igg.rds")
type <- "IGG"

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


# do a pca to check if the correct plate is loaded

pca <- prcomp(x = plate[,2:ncol(plate)],scale. = T,center = T)

percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                       condition = as.factor(labels))

g2 <- ggplot(data_pca, aes(PC1, PC2)) +
  geom_point(aes(colour = condition))  +
  ggtitle("PCA plot of the samples") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("darkorange2", "dodgerblue4","green","red","yellow"))
g2


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

#############################################################

# calculate negative control threshold
cutoff <- calc_thresh(ncov,outlier = 1.5)


# roc analysis
# set weights to sensitivity for Youden point calculation
w <- c(0.25,0.25,0.25)

#set output dir
output_dir <- "performance_report_aggregated/igg/ratios_7ug"

roc_perf <- roc_anal(ncov = ncov,n = 1000,output_dir = paste0(output_dir,"/roc_analysis"),conc = conc,type = type,w=w)

# add barplot and scatterplot
plot_scatter_bar(ncov = ncov_all,roc_perf = roc_perf,cutoff = cutoff,
                 output_dir = paste0(output_dir,"/plots"),
                 type = type,conc = conc)

# evaluate performance roc thresh
performance_roc <- eval_performance_roc(ncov = ncov,roc_perf = roc_perf,
                     type = type,
                     conc = conc[3],
                     output_dir = paste0(output_dir,"/performance_eval"),
                     normal = "raw")

# evaluate performance negative control thresh
performance_thresh <- eval_performance_thresh(ncov = ncov,cutoff = cutoff,roc_perf = roc_perf,
                        type = type,
                        conc = conc[3],
                        output_dir = paste0(output_dir,"/performance_eval"),
                        normal = "raw")

# write overall performance with CI
performance <- bind_rows(performance_roc,performance_thresh)

#write results
write.csv(performance,paste0(output_dir,"/performance_report.csv"))
