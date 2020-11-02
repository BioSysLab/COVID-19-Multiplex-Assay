
abbot <- read.csv("abbott_epidem/abbot.csv")

match <- read.csv("abbott_epidem/match.csv",header = F)
colnames(match) <- c("abbott","prot")

abbot_cov <- abbot[grepl(pattern = "COVID",x = abbot$sample),]
abbot_nocov <- abbot[!grepl(pattern = "COVID",x = abbot$sample),]

length(which(abbot_cov$sample %in% match$abbott))

abbot_cov <- abbot_cov[which(abbot_cov$sample %in% match$abbott),]

abbot_cov <- left_join(abbot_cov,match,by = c("sample"="abbott"))

abbot_cov <- abbot_cov[,c(3,2)]

colnames(abbot_cov)<-c("sample","value")

abbot_new <- bind_rows(abbot_nocov,abbot_cov)

length(which(abbot_new$sample %in% ncov$Sample))

abbot_new <- abbot_new[which(abbot_new$sample %in% ncov$Sample),]

new_ncov <- left_join(abbot_new,ncov,by=c("sample"="Sample"))

ab_plot <- ggplot(new_ncov,aes(x=(N),y = (value))) +
  geom_point(size = 1.2) + xlab("N multiplex") + ylab("N Abbott") + 
  scale_x_continuous(breaks = c(0,5.48,10,15),labels=c("0" = "0", "5.48" = "Cut-off = 5.48","10" = "10","15"="15"))+
  scale_y_continuous(breaks = c(0,1.4,2.5,5,7.5),labels=c("0" = "0.0", "1.4" = "Cut-off = 1.4","2.5" = "2.5","5" = "5.0","7.5"="7.5"))+
  geom_vline(aes(xintercept=5.48),linetype="dotdash")+
  geom_hline(aes(yintercept=1.4),linetype="dotdash")+
  annotate(geom = 'text', label = 'Pearson r = 0.20', x = 3, y = 3, hjust = 0, vjust = 1,size = 7)+
  ggtitle("Total (IgG/IgA/IgM) N Antigen") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal(base_size = 20)
png(filename = paste0("figures","/","epidem/","abbott_scatterplot.png"),width = 8,height = 6,units = "in",res = 600)
ab_plot
dev.off()

ab_plot <- ggplot(new_ncov,aes(x=(N),y = (value))) +
  geom_point(size = 1.2) + xlab("N multiplex") + ylab("N Abbott") +ylim(c(0,7.5)) + xlim(c(0,15))+
  geom_vline(aes(xintercept=5.48),linetype="dotdash")+
  geom_hline(aes(yintercept=1.4),linetype="dotdash")+
  annotate(geom = 'text', label = 'Pearson r = 0.20', x = 3, y = 3, hjust = 0, vjust = 1,size = 7)+
  ggtitle("Total (IgG/IgA/IgM) N Antigen") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal(base_size = 20)

ab_plot+scale_y_continuous(breaks = sort(c(ab_plot$layout$coord$labels(ab_plot$layout$panel_params)[[1]]$y.major_source,1.4)))
