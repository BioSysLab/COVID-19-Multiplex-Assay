library(tidyverse)
library(viridis)
library(reshape2)
plate1 <- read.csv("performance_b/plate1_igm.csv")
#plate1 <- plate1[,-1]
plate1$Sample <- as.character(plate1$Sample)

plate2 <- read.csv("performance_b/plate2_igm.csv")
plate2$Sample <- as.character(plate2$Sample)

labels1 <- c(rep("Negative",78),rep("Positive",15))
labels2 <- c(rep("Blind samples",40),rep("Positive",53))
labels <- c(labels1,labels2)
labels2 <- c(rep("negative_gennimatas",66),rep("negative_patras",12),rep("positive_patras",15),rep("test",40),rep("positive_alexandras",53))


igg <- bind_rows(plate1,plate2)
igg$labels <- labels

correlations <- cor(igg[,2:20])

rownames(correlations) <- colnames(igg)[2:20]
colnames(correlations) <- colnames(igg)[2:20]

# heatmap of feature correlations

melted_cor <- melt(correlations)

plot_cors <- ggplot(data =melted_cor, aes(x=Var1, y=Var2, fill= value)) + 
  geom_tile()+scale_fill_viridis(limits = c(-1,1),direction = 1) + 
  labs( fill = "correlation") + xlab("antigen")+ylab("antigen")+
  ggtitle("antigen correlation") +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

png(file=paste0("iga_results/correlation_plot_iga.png"),width=8,height=6,units = "in",res=300)
print(plot_cors)
dev.off()

# pca 

pca <- prcomp(x = igg[,2:20],scale. = T,center = T)

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
# write pca plot
png(file=paste0("iga_results/pca_plot_iga_allvars.png"),width=8,height=6,units = "in",res=300)
print(g2)
dev.off()
#
#pca of old hcov

hcov <- c("X40021.V08H_30UG","X40600.V08H_30UG","X40601.V08H_30UG")
old_cov <- igg[,c(1,21,which(colnames(igg) %in% hcov))]

pca_old <- prcomp(x = old_cov[,3:5],scale. = T,center = T)

percentVar <- round(100*pca_old$sdev^2/sum(pca_old$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_pca_old <- data.frame(PC1 = pca_old$x[,1], PC2 = pca_old$x[,2],
                       condition = as.factor(labels))

g2_old <- ggplot(data_pca_old, aes(PC1, PC2)) +
  geom_point(aes(colour = condition))  +
  ggtitle("PCA plot of the samples for old cov") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("darkorange2", "dodgerblue4","green","red","yellow"))
g2_old

# write pca plot
png(file=paste0("iga_results/pca_plot_iga_oldcov.png"),width=8,height=6,units = "in",res=300)
print(g2_old)
dev.off()

# pca of new cov features

ncov <- c("REC31828.100_30UG","X40592.V08H_30UG","REC31812.100_30UG","REC31806.500_30UG")
ncov <- igg[,c(1,21,which(colnames(igg) %in% ncov))]
colnames(ncov) <- c("Sample","labels","Nucleoprotein","Spike glycoprotein (1st)","Spike RBD protein","Spike glycoprotein (2nd)")
pca_new <- prcomp(x = ncov[,3:ncol(ncov)],scale. = T,center = T)

percentVar <- round(100*pca_new$sdev^2/sum(pca_new$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_pca_new <- data.frame(PC1 = pca_new$x[,1], PC2 = pca_new$x[,2],
                           condition = as.factor(labels))

g2_new <- ggplot(data_pca_new, aes(PC1, PC2)) +
  geom_point(aes(colour = condition))  +
  ggtitle("PCA plot of the samples for new cov") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("darkorange2", "dodgerblue4","green","red","yellow"))
g2_new

# write pca plot
png(file=paste0("iga_results/pca_plot_iga_newcov.png"),width=8,height=6,units = "in",res=300)
print(g2_new)
dev.off()


#write general figures

igg$labels <- labels
######## outlier parameter #############
outlier <- 1.5
for (i in 2:20) {
  meas <- as.numeric(igg[,i])
  thresh1 <- mean(meas[labels=="Negative"]) + 3* sd(meas[labels=="Negative"])
  id <- which(meas[labels=="Negative"]>(outlier)*thresh1)
  if (length(id) != 0) {
    neg <- which(labels=="Negative")
    neg <- neg[-id]
    thresh1 <- mean(meas[neg]) + 3* sd(meas[neg])
  }
  jpeg(file=paste0("iga_results/iga_",colnames(igg)[i],".jpeg"),width=8,height=6,units = "in",quality=150,res=300)
  general <- ggplot(igg,aes(1:186,igg[,i]))+geom_point(aes(colour = labels))+xlab("Samples")+ylab("Fluorescent Intensity") + geom_hline(aes(yintercept=thresh1, linetype="cut-off"), color = "red")+
    scale_linetype_manual(name = "cutoff", values = 2, 
                                  guide = guide_legend(override.aes = list(color = "red")))+
    annotate(geom="text", x=10, y=(-5000), label="Negatives",color="black")+
    annotate(geom="text", x=10, y=(thresh1+10000), label="Positives",color="black")+ggtitle(colnames(igg)[i])
  gen <- ggplot_build(general)
  general <- general + scale_y_continuous(breaks = sort(c(gen$layout$coord$labels(gen$layout$panel_params)[[1]]$y.major_source,thresh1)))
  print(general)
  dev.off()
  
  #df <- igg[,c(1,i,21)]
  #barplot <- ggplot(df,aes(x=reorder(Sample,-df[,2]), y = df[,2],width=.5,fill = labels)) +    
   # geom_bar(stat="identity") +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5))+
  #  xlab("Samples") + ylab("Fluorescent Intensity") + ggtitle(colnames(igg)[i])+ 
   # geom_hline(aes(yintercept = thresh1, linetype="cut-off"), color = "red")
  #bar <- ggplot_build(barplot)
  #barplot <- barplot + scale_y_continuous(breaks = sort(c(bar$layout$coord$labels(bar$layout$panel_params)[[1]]$y.major_source,thresh1)))
  #jpeg(file=paste0("igg_results/igg_",colnames(igg)[i],".jpeg"),width=8,height=6,units = "in",quality=150,res=300)
  #print(barplot)
}

#write new cov scatter plots
excl <- c("48501","48524","48530")
ncov <- ncov[-which(ncov$Sample%in%excl),]
gscatters <- NULL
k <- 1
for (i in 3:ncol(ncov)) {
  meas <- as.numeric(ncov[,i])
  thresh1 <- mean(meas[ncov$labels=="Negative"]) + 3* sd(meas[ncov$labels=="Negative"])
  id <- which(meas[ncov$labels=="Negative"]>(outlier)*thresh1)
  if (length(id) != 0) {
    neg <- which(ncov$labels=="Negative")
    neg <- neg[-id]
    thresh1 <- mean(meas[neg]) + 3* sd(meas[neg])
  }
  #jpeg(file=paste0("igg_results/igg_",colnames(igg)[i],".jpeg"),width=8,height=6,units = "in",quality=150,res=300)
  general <- ggplot(ncov,aes(1:183,meas))+geom_point(aes(colour = labels))+xlab("Samples")+ylab("Fluorescent Intensity") + geom_hline(aes(yintercept=thresh1, linetype="cut-off"), color = "red")+
    scale_linetype_manual(name = "cutoff", values = 2, 
                          guide = guide_legend(override.aes = list(color = "red")))+
    annotate(geom="text", x=10, y=(-5000), label="Negatives",color="black")+
    annotate(geom="text", x=10, y=(thresh1+10000), label="Positives",color="black")+ggtitle(colnames(ncov)[i])
  gen <- ggplot_build(general)
  general <- general + scale_y_continuous(breaks = sort(c(gen$layout$coord$labels(gen$layout$panel_params)[[1]]$y.major_source,thresh1)))
  gscatters[[k]] <- ggplotGrob(general)
  general <- NULL
  gen <- NULL
  k <- k+1
  #print(general)
  #dev.off()
}
jpeg(file=paste0("iga_results/ncov_antigen_scatterplot.jpeg"),width=12,height=9,units = "in",res=300)
gridExtra::grid.arrange(grobs=gscatters,nrow=2,top="Antigen Performance Analysis")
dev.off()

#write barplot
for (i in 2:20) {
  meas <- as.numeric(igg[,i])
  thresh1 <- mean(meas[labels=="Negative"]) + 3* sd(meas[labels=="Negative"])
  id <- which(meas[labels=="Negative"]>(outlier)*thresh1)
  if (length(id) != 0) {
    neg <- which(labels=="Negative")
    neg <- neg[-id]
    thresh1 <- mean(meas[neg]) + 3* sd(meas[neg])
  }
  df <- igg[,c(1,i,21)]
  barplot <- ggplot(df,aes(x=reorder(Sample,-df[,2]), y = df[,2],width=0.5,fill = labels)) +    
    geom_bar(stat="identity",position=position_nudge((x=0.5))) +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=1,size=2))+
    xlab("Samples") + ylab("Fluorescent Intensity") + ggtitle(colnames(igg)[i])+ 
    geom_hline(aes(yintercept = thresh1, linetype="cut-off"), color = "red")
  bar <- ggplot_build(barplot)
  barplot <- barplot + scale_y_continuous(breaks = sort(c(bar$layout$coord$labels(bar$layout$panel_params)[[1]]$y.major_source,thresh1)))
  jpeg(file=paste0("iga_results/iga_bar_",colnames(igg)[i],".jpeg"),width=9,height=6,units = "in",quality=150,res=300)
  print(barplot)
  dev.off()
}
# write ncov barplots
gbars <- NULL
k <- 1
for (i in 3:ncol(ncov)) {
  meas <- as.numeric(ncov[,i])
  thresh1 <- mean(meas[ncov$labels=="Negative"]) + 3* sd(meas[ncov$labels=="Negative"])
  id <- which(meas[ncov$labels=="Negative"]>(outlier)*thresh1)
  if (length(id) != 0) {
    neg <- which(ncov$labels=="Negative")
    neg <- neg[-id]
    thresh1 <- mean(meas[neg]) + 3* sd(meas[neg])
  }
  df <- ncov[,c(1,i,2)]
  barplot <- ggplot(df,aes(x=reorder(Sample,-df[,2]), y = df[,2],width=0.5,fill = labels)) +    
    geom_bar(stat="identity",position=position_nudge((x=0.5))) +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=1,size=2))+
    xlab("Samples") + ylab("Fluorescent Intensity") + ggtitle(colnames(ncov)[i])+ 
    geom_hline(aes(yintercept = thresh1, linetype="cut-off"), color = "red")
  bar <- ggplot_build(barplot)
  barplot <- barplot + scale_y_continuous(breaks = sort(c(bar$layout$coord$labels(bar$layout$panel_params)[[1]]$y.major_source,thresh1)))
  #jpeg(file=paste0("igg_iga_igm_results/ig_all_bar_",colnames(igg)[i],".jpeg"),width=9,height=6,units = "in",quality=150,res=300)
  #print(barplot)
  #dev.off()
  gbars[[k]] <- ggplotGrob(barplot)
  k <- k+1
  barplot <- NULL
  bar <- NULL
}

jpeg(file=paste0("iga_results/ncov_antigen_barplot.jpeg"),width=12,height=9,units = "in",res=300)
gridExtra::grid.arrange(grobs=gbars,nrow=2,top="Antigen Performance Analysis")
dev.off()

# sensitivity, specificity, accuracy, ppv, npv, roc auc

library(ROCit)
library(tidyverse)
igg_val <- igg %>% filter(labels != "Blind samples")
#meas <- as.numeric(igg_val$REC31828.100_15UG)
#labels <- igg_val$labels
evaluate_antigen <- function(meas,labels,antigen_name,output_dir,outlier){
  library(ROCit)
  roc_empirical <- rocit(score = meas, class = as.factor(labels),
                         negref = "Negative")
  auc <- ciAUC(roc_empirical)
  roc_ci <- ciROCcustom(roc_empirical, level = 0.95)
  png(file=paste0(output_dir,"/roc_curve","_",antigen_name,".png"),width=8,height=6,units = "in",res=300)
  plot(roc_empirical,values=F)
  lines(roc_ci$FPR,roc_ci$LowerTPR,col= "red",lty =2)
  lines(roc_ci$FPR,roc_ci$UpperTPR,col= "red",lty =2)
  dev.off()
  save_auc <- c(auc$lower,auc$AUC,auc$upper,auc$conf_level)
  names(save_auc) <- c("lower bound","AUC","upper bound","conf_level")
  #write.csv(save_auc,paste0("auc_ci_",antigen_name,".csv"))
  roc_anal <- cbind(roc_empirical$Cutoff,roc_empirical$TPR,1-roc_empirical$FPR)
  colnames(roc_anal) <- c("threshold","sensitivity","specificity")
  write.csv(roc_anal,paste0(output_dir,"/roc_analysis_",antigen_name,".csv"))
  perf_metrics <- function(meas,labels,thresh){
    df <- data.frame(meas,labels)
    colnames(df) <- c("meas","labels")
    tp <- df %>% filter(meas>=thresh) %>% filter(labels=="Positive")
    tp <- nrow(tp)
    
    fp <- df %>% filter(meas>=thresh) %>% filter(labels=="Negative")
    fp <- nrow(fp)
    
    tn <- df %>% filter(meas<=thresh) %>% filter(labels=="Negative")
    tn <- nrow(tn)
    
    fn <- df %>% filter(meas<=thresh) %>% filter(labels=="Positive")
    fn <- nrow(fn)
    
    sens <- tp/(tp+fn)
    spec <- tn/(tn+fp)
    ppv <- tp/(tp+fp)
    npv <- tn/(tn+fn)
    acc <- (tn+tp)/(tn+fn+tp+fp)
    # add confidence intervals
    return(c(sens, spec, ppv,npv, acc))
  }
  
  # calculate threshold alexopouleiros
  thresh1 <- mean(meas[labels=="Negative"]) + 3* sd(meas[labels=="Negative"])
  id <- which(meas[labels=="Negative"]>(outlier)*thresh1)
  if (length(id) != 0) {
    neg <- which(labels=="Negative")
    neg <- neg[-id]
    thresh1 <- mean(meas[neg]) + 3* sd(meas[neg])
  }
  metrics1 <- perf_metrics(meas,labels,thresh1)
  #write.csv(metrics1,paste0("performance_metrics_",antigen_name,".csv"))
  #threshold nikopoulos
  box <- boxplot(meas[labels=="Positive"])
  thresh2 <- box$stats[2,1]
  metrics2 <- perf_metrics(meas,labels,thresh2)
  return(c(metrics1,save_auc))
}

results <- NULL
for (i in 2:20) {
  metrics <- evaluate_antigen(meas = as.numeric(igg_val[,i]),labels = igg_val$labels,antigen_name = colnames(igg_val)[i],
                              output_dir = "iga_results",outlier=outlier)
  results <- rbind(results,metrics)
}
colnames(results) <- c("sensitivity","specificity","ppv","npv","accuracy","lower bound","AUC","upper bound","conf_level")
rownames(results) <- colnames(igg_val)[2:20]
results <- as.data.frame(results)
write.csv(results,"iga_results/performance_metrics_iga.csv")



library(caret)
library(tidyverse)
#decision tree
general_cov <- c("REC31828.100_7.5UG","X40592.V08H_30UG","REC31806.500_30UG","X40021.V08H_30UG","X40600.V08H_30UG","X40601.V08H_30UG")
gen_cov <- igg[,c(1,21,which(colnames(igg) %in% general_cov))]

excl <- c("")

dataset <- ncov %>% column_to_rownames("Sample") %>% filter(labels != "Blind samples")
#dataset[,2:ncol(dataset)] <- scale(dataset[,2:ncol(dataset)])

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

alexo_tree = train(labels ~ ., 
                  data=dataset, 
                  method="rpart", 
                  trControl = fitControl,metric = "ROC",tuneLength = 10)
alexo_tree

par(xpd = NA) # Avoid clipping the text in some device
plot(alexo_tree$finalModel)
text(alexo_tree$finalModel,  digits = 3)
