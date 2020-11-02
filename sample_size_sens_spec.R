library(binomSamSize)
library(tidyverse)


wilson <- function(p,CL,n){
  if (CL == 0.95) {
    z <- 1.96
  }
  if (CL == 0.99) {
    z <- 2.576
  }
  if (CL == 0.90) {
    z <- 1.645
  }
  p_dot <- (p*n+(z^2)/2)/(n+z^2)
  n_dot <- n+z^2
  interval_half <- (z/n_dot)*sqrt(n*p*(1-p)+z^2/4)
  return(list(p_dot,interval_half))
}

sens <- seq(0.9,0.9,0.05)
n_pos <- seq(135,150,1)
results_sensitivity <- matrix(0,nrow = length(n_pos)*length(sens),ncol = 3)
results_sensitivity <- as.data.frame(results_sensitivity)
colnames(results_sensitivity) <- c("positives","anticipated_sensitivity","half_width_95")
row <- 1
for (i in 1:length(n_pos)) {
    for (k in 1:length(sens)) {
      scores <- wilson(sens[k],CL = 0.95,n = n_pos[i])
      results_sensitivity$half_width_95[row] <- scores[[2]]
      results_sensitivity$positives[row] <- n_pos[i]
      results_sensitivity$anticipated_sensitivity[row] <- sens[k]
      row <- row+1
    }
}
results_sensitivity$anticipated_sensitivity <- as.factor(results_sensitivity$anticipated_sensitivity)



png(file="result_sens95.png",width=7,height=6,units = "in",res=300)
plot(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.70],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity==0.70],
     xlim = c(9,41),ylim = c(0,0.25),"o",col = "black",lwd = 1.7,xaxs = "i",lty = 5,pch = 1,yaxs="i",cex = 0.6,
     xlab = "number of positive cases",ylab = "95% interval half width")
lines(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.75],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity==0.75],
      col = "#e41a1c",lwd = 1.7,type = "o",lty = 5,pch = 2,cex = 0.6)
lines(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.80],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity==0.80],
      col = "yellow",lwd = 1.7,type = "o",lty = 5,pch = 3,cex = 0.6)
lines(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.85],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity==0.85],
      col = "#377eb8",lwd = 1.7, type = "o",lty = 5,pch = 4,cex = 0.6)
lines(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.90],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity== 0.90],
      col = "green",lwd = 1.7, type = "o",lty = 5 , pch = 5,cex = 0.6)
lines(results_sensitivity$positives[results_sensitivity$anticipated_sensitivity==0.95],results_sensitivity$half_width_95[results_sensitivity$anticipated_sensitivity== 0.95],
      col = "purple",lwd = 1.7, type = "o",lty = 5 , pch = 6,cex = 0.6)
abline(v=25,col="red",lwd = 1)
title("Estimated \u00b1 intervals (95% confidence) for sensitivity",adj = 0)
legend("topright", 
       legend = c("sensitivity = 0.70", "sensitivity = 0.75", "sensitivity = 0.80","sensitivity = 0.85","sensitivity = 0.90","sensitivity = 0.95"),
       col = c('black', "#e41a1c","yellow",
               '#377eb8',"green","purple"),pch =c(1,2,3,4,5,6),
       pt.cex = 0.6, lwd = 0.8,
       cex = 0.6, 
       text.col = "black", 
       horiz = F )
dev.off()

spec <- seq(0.85,0.95,0.05)
n_neg <- seq(40,100,10)
results_specificity <- matrix(0,nrow = length(n_neg)*length(spec),ncol = 3)
results_specificity <- as.data.frame(results_specificity)
colnames(results_specificity) <- c("negatives","anticipated_specificity","half_width_95")
row <- 1
for (i in 1:length(n_neg)) {
  for (k in 1:length(spec)) {
    scores <- wilson(spec[k],CL = 0.95,n = n_neg[i])
    results_specificity$half_width_95[row] <- scores[[2]]
    results_specificity$negatives[row] <- n_neg[i]
    results_specificity$anticipated_specificity[row] <- spec[k]
    row <- row+1
  }
}
results_specificity$anticipated_specificity <- as.factor(results_specificity$anticipated_specificity)

t1 <- results_specificity %>% filter(negatives==40)
t2 <- results_specificity %>% filter(negatives==80)

res <- 100* (t1$half_width_95 - t2$half_width_95 )/t1$half_width_95

png(file="result_spec90.png",width=7,height=6,units = "in",res=300)
plot(results_specificity$negatives[results_specificity$anticipated_specificity==0.70],results_specificity$half_width_95[results_specificity$anticipated_specificity==0.70],
     xlim = c(39,101),ylim = c(0,0.15),"o",col = "black",lwd = 1.7,xaxs = "i",lty = 5,pch = 1,yaxs="i",cex = 0.6,
     xlab = "number of negative cases",ylab = "95% interval half width")
lines(results_specificity$negatives[results_specificity$anticipated_specificity==0.75],results_specificity$half_width_95[results_specificity$anticipated_specificity==0.75],
      col = "#e41a1c",lwd = 1.7,type = "o",lty = 5,pch = 2,cex = 0.6)
lines(results_specificity$negatives[results_specificity$anticipated_specificity==0.80],results_specificity$half_width_95[results_specificity$anticipated_specificity==0.80],
      col = "yellow",lwd = 1.7,type = "o",lty = 5,pch = 3,cex = 0.6)
lines(results_specificity$negatives[results_specificity$anticipated_specificity==0.85],results_specificity$half_width_95[results_specificity$anticipated_specificity==0.85],
      col = "#377eb8",lwd = 1.7, type = "o",lty = 5,pch = 4,cex = 0.6)
lines(results_specificity$negatives[results_specificity$anticipated_specificity==0.90],results_specificity$half_width_95[results_specificity$anticipated_specificity== 0.90],
      col = "green",lwd = 1.7, type = "o",lty = 5 , pch = 5,cex = 0.6)
lines(results_specificity$negatives[results_specificity$anticipated_specificity==0.95],results_specificity$half_width_95[results_specificity$anticipated_specificity== 0.95],
      col = "purple",lwd = 1.7, type = "o",lty = 5 , pch = 6,cex = 0.6)
abline(v=80,col="red",lwd = 1)
title("Estimated \u00b1 intervals (90% confidence) for specificity",adj = 0)

legend("topright", 
       legend = c("specificity = 0.70", "specificity = 0.75", "specificity = 0.80","specificity = 0.85","specificity = 0.90","specificity = 0.95"),
       col = c('black', "#e41a1c","yellow",
               '#377eb8',"green","purple"),pch =c(1,2,3,4,5,6),
       pt.cex = 0.6, lwd = 0.8,
       cex = 0.6, 
       text.col = "black", 
       horiz = F )
dev.off()



library(pwr)


pwr.p.test(h = ES.h(p1 = 0.95, p2 = 0.90),
           n=20,
           sig.level = 0.05,
           power = NULL,
           alternative = "greater")


pwr.2p.test(h = ES.h(p1 = 0.90, p2 = 0.75), sig.level = 0.05, power = NULL, n = 100)


install.packages("DTComPair")



library(DTComPair)




num <- 1
dn <- rep(0,80)
n <- c(10,20,30,40,50)
sens <- c(0.7,0.8,0.9)
tab_sens <- matrix(0,nrow = length(n)*length(sens)^2,ncol = 4)
tab_sens <- as.data.frame(tab_sens)
colnames(tab_sens) <- c("positives","sensitivity_1","sensitivity_2","p_value")
for (k in 1:length(n)) {
  for (i in 1:length(sens)) {
    for (j in 1:length(sens)) {
      dp <- rep(1,n[k])
      d <- c(dp,dn)
      y1n <- rep(0,80)
      y2n <- rep(0,80)
      y1pp <- rep(1,n[k]*sens[i])
      y1pn <- rep(0,n[k]-length(y1pp))
      y1p <- c(y1pp,y1pn)
      
      y2pp <- rep(1,n[k]*sens[j])
      y2pn <- rep(0,n[k]-length(y2pp))
      y2p <- c(y2pp,y2pn)
      
      y1 <- c(y1p,y1n)
      y2 <- c(y2p,y2n)
      
      tab <- tab.paired(d, y1, y2, data = NULL, c("test1","test2"))
      test <- sesp.mcnemar(tab)
      tab_sens$positives[num] <- n[k]
      tab_sens$sensitivity_1[num] <- test$sensitivity$test1
      tab_sens$sensitivity_2[num] <- test$sensitivity$test2
      tab_sens$p_value[num] <- test$sensitivity$p.value
      num <- num+1
    }
  }
}




num <- 1
dp <- rep(1,20)
n <- c(40,50,60,70,80,90,100)
spec <- c(0.7,0.8,0.9)
tab_spec <- matrix(0,nrow = length(n)*length(spec)^2,ncol = 4)
tab_spec <- as.data.frame(tab_spec)
colnames(tab_spec) <- c("negatives","specificity_1","specificity_2","p_value")
for (k in 1:length(n)) {
  for (i in 1:length(spec)) {
    for (j in 1:length(spec)) {
      dn <- rep(0,n[k])
      d <- c(dp,dn)
      y1p <- rep(1,20)
      y2p <- rep(1,20)
      y1nn <- rep(0,n[k]*spec[i])
      y1np <- rep(1,n[k]-length(y1nn))
      y1n <- c(y1nn,y1np)
      
      y2nn <- rep(0,n[k]*spec[j])
      y2np <- rep(1,n[k]-length(y2nn))
      y2n <- c(y2nn,y2np)
      
      y1 <- c(y1p,y1n)
      y2 <- c(y2p,y2n)
      
      tab <- tab.paired(d, y1, y2, data = NULL, c("test1","test2"))
      test <- sesp.mcnemar(tab)
      tab_spec$negatives[num] <- n[k]
      tab_spec$specificity_1[num] <- test$specificity$test1
      tab_spec$specificity_2[num] <- test$specificity$test2
      tab_spec$p_value[num] <- test$specificity$p.value
      num <- num+1
    }
  }
}

tab_sens <- tab_sens %>% filter(!is.na(p_value))
tab_sens <- tab_sens %>% filter(sensitivity_1 > sensitivity_2)
tab_spec <- tab_spec %>% filter(!is.na(p_value))
tab_spec <- tab_spec %>% filter(specificity_1 > specificity_2)

write.csv(tab_sens,"paired_testing_between_diagnostic_tests_sample_size_sensitivity.csv")
write.csv(tab_spec,"paired_testing_between_diagnostic_tests_sample_size_specificity.csv")
