library(tidyverse)
library(viridis)
library(reshape2)
library(doFuture)
library(ggpubr)
library(ggrepel)

performance_total <- read.csv("paper_results/total/total4_wilson_performance.csv")
performance_igg <- read.csv("paper_results/igg/igg4_wilson_performance.csv")
performance_iga <- read.csv("paper_results/iga/iga4_wilson_performance.csv")
performance_igm <- read.csv("paper_results/igm/igm4_wilson_performance.csv")

performance <- bind_rows(performance_total,performance_igg,performance_iga,performance_igm)

# fix names

performance$type[which(performance$type == "total")] <- "Total (IgG/IgA/IgM)"
performance$type[which(performance$type == "igg")] <- "IgG"
performance$type[which(performance$type == "iga")] <- "IgA"
performance$type[which(performance$type == "igm")] <- "IgM"

performance$rule <- as.character(performance$rule)
performance$rule[which(performance$rule == "RBD only")] <- "RBD"
performance$rule[which(performance$rule == "S1 only")] <- "S1"
performance$rule[which(performance$rule == "N only")] <- "N"




# fix upper thresholds

performance$sensitivity_upper[performance$sensitivity_upper>1] <- 1
performance$specificity_upper[performance$specificity_upper>1] <- 1
performance$accuracy_upper[performance$accuracy_upper>1] <- 1

#colnames(performance)[which(colnames(performance)=="rule")] <- "SARS-CoV-2 antigen"

#performance_single <- performance %>% filter(rule == "RBD" | rule == "S1" |rule=="N")

# plot singles
library("ggsci")
library("viridis")
performance_plot <- performance %>% 
  filter(rule == "RBD" | rule == "S1" |rule=="N"| rule == "rule" | rule == "all" | rule == "Majority")
performance_plot$rule[which(performance_plot$rule == "all")] <- "RBD + N + S1"
performance_plot$rule[which(performance_plot$rule == "rule")] <- "RBD + N/S1"

#fix order of type

performance_plot$type <- factor(performance_plot$type, levels = c("IgM","IgA","IgG","Total (IgG/IgA/IgM)"))

specificity_single <- ggplot(data = performance_plot, aes(x = type, y = specificity, fill = rule))+
  geom_bar(stat="identity", position=position_dodge(0.9))+
  geom_text(aes(y=rep(0.925, times = 24), label=rule), 
  size = 5,color = "white",
  position = position_dodge(width=0.9) 
  ) +
  theme_minimal(base_size = 12)+xlab("")+ylab("Specificity")+labs(fill = "SARS-CoV-2 antigen")+
  geom_errorbar(aes(ymin=specificity_lower, ymax=specificity_upper), width=.2,
                position=position_dodge(.9))+
  coord_flip(ylim = c(0.9,1))+theme(axis.text.y = element_text(angle = 90,hjust = 0.5,face = "bold",size = 12),
                                    axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
                                    legend.text = element_text(size = 12))+
  scale_fill_manual(values = c("#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02"))
specificity_single  

sensitivity_single <- ggplot(data = performance_plot, aes(x = type, y = sensitivity, fill = rule))+
  geom_bar(stat="identity", position=position_dodge(0.9))+
  geom_text(aes(y=rep(0.3, times = 24), label=rule), 
            size = 5,color = "white",
            position = position_dodge(width=0.9) 
  ) +
  theme_minimal(base_size = 12)+xlab("")+ylab("Sensitivity")+labs(fill = "SARS-CoV-2 antigen")+
  geom_errorbar(aes(ymin=sensitivity_lower, ymax=sensitivity_upper), width=.2,
                position=position_dodge(.9))+
  coord_flip(ylim = c(0.1,1))+theme(axis.text.y = element_text(angle = 90,hjust = 0.5,face = "bold",size = 12),
                                    axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
                                    legend.text = element_text(size = 12))+
  scale_fill_manual(values = c("#1b9e77",
                               "#d95f02",
                               "#7570b3",
                               "#e7298a",
                               "#66a61e",
                               "#e6ab02"))
sensitivity_single  
