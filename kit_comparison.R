library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
library(ggpubr)

abbot <- read.csv("data_exp2/abbot.csv")
euro <- read.csv("data_exp2/euroimmune.csv")

ab_plot <- ggplot(abbot,aes(x=(igg_N),y = (abbott))) +
  geom_point(size = 0.5) + xlab("N multiplex") + ylab("N Abbott") + 
  scale_x_continuous(breaks = c(0,3.24,10,15),labels=c("0" = "0", "3.24" = "Cut-off = 3.24","10" = "10","15"="15"))+
  scale_y_continuous(breaks = c(0,1.4,2.5,5,7.5),labels=c("0" = "0.0", "1.4" = "Cut-off = 1.4","2.5" = "2.5","5" = "5.0","7.5"="7.5"))+
  geom_vline(aes(xintercept=3.24),linetype="dotdash")+
  geom_hline(aes(yintercept=1.4),linetype="dotdash")+
  annotate(geom = 'text', label = 'Pearson r = 0.98', x = 10, y = 7.5, hjust = 0, vjust = 1,size = 7)+
  ggtitle("IgG N Antigen") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal(base_size = 20)
png(filename = paste0("figures","/","validation/","abbott_scatterplot.png"),width = 8,height = 6,units = "in",res = 600)
ab_plot
dev.off()


euro_plot <- ggplot(euro,aes(x=(igg_s1),y = (euro))) +
  geom_point(size = 3) + xlab("S1 multiplex") + ylab("S1 Euroimmun") + 
  scale_x_continuous(breaks = c(0,5.22,10,15,20),labels=c("0" = "0", "5.22" = "Cut-off = 5.22","10" = "10","15"="15","20"="20"))+
  scale_y_continuous(breaks = c(0,1.1,5,10,15,20,25),labels=c("0" = "0", "1.1" = "Cut-off = 1.1","5" = "5","10" = "10","15"="15",
                                                              "20"="20","25"="25"))+
  geom_vline(aes(xintercept=5.22),linetype="dotdash")+
  geom_hline(aes(yintercept=1.1),linetype="dotdash")+
  annotate(geom = 'text', label = 'Pearson r = 0.90', x = 10, y = 15, hjust = 0, vjust = 1,size = 7)+
  ggtitle("IgG S1 Antigen") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal(base_size = 20)
png(filename = paste0("figures","/","validation/","euroimmune_scatterplot.png"),width = 8,height = 6,units = "in",res = 600)
euro_plot
dev.off()

png(filename = paste0("figures","/","validation/","scatterplot.png"),width = 17,height = 7,units = "in",res = 600)
combined <- ggarrange(euro_plot,ab_plot, ncol=2, nrow=1)
annotate_figure(combined)
dev.off()
