library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
library(ggpubr)

euro <- read.csv("validation/euroimmune.csv")
colnames(euro) <- c("S1 Euroimmun","S1 Multiplex","N Multiplex","RBD Multiplex","RBD & N|S1 rule")
euro$idx <- 1:60
euro_melted <- melt(euro,id.vars = "idx")
euro_melted$variable <- factor(euro_melted$variable, levels = c("RBD Multiplex","N Multiplex","S1 Multiplex","RBD & N|S1 rule","S1 Euroimmun"))
euro_melted$value <- as.factor(euro_melted$value)


euro_plot <- ggplot(data = euro_melted, aes(x=variable, y=idx)) + 
  geom_tile(aes(fill = value,width=0.8, height=0.8),show.legend = T) +
  scale_fill_discrete(name="Result",
                      labels=c("Negative","Positive"))+xlab("")+
  xlab("")+
  scale_x_discrete(position = "top")+
  ylab("Samples") + theme(axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(),text = element_text(size = 15))


abbott <- read.csv("validation/abbot.csv")
abbott <- sample_n(tbl = abbott,size = 31)
colnames(abbott) <- c("N Abbott","N Multiplex","S1 Multiplex","RBD Multiplex","RBD & N|S1 rule")
abbott$idx <- 1:31
abbott_melted <- melt(abbott,id.vars = "idx")
abbott_melted$variable <- factor(abbott_melted$variable, levels = c("RBD Multiplex","N Multiplex","S1 Multiplex","RBD & N|S1 rule","N Abbott"))
abbott_melted$value <- as.factor(abbott_melted$value)

abbott_plot <- ggplot(data = abbott_melted, aes(x=variable, y=idx)) + 
  geom_tile(aes(fill = value,width=0.8, height=0.8),show.legend = T) +
  scale_fill_discrete(name="Result",
                       labels=c("Negative","Positive"))+xlab("")+
  scale_x_discrete(position = "top")+
  ylab("Samples") + theme(axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(),text = element_text(size = 15))


png(filename = paste0("figures","/","validation/","heatmap_legend.png"),width = 17,height = 7,units = "in",res = 600)
combined <- ggarrange(euro_plot,abbott_plot, ncol=2, nrow=1,common.legend = TRUE, legend="bottom")
annotate_figure(combined)
dev.off()
