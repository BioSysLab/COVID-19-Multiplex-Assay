library(binomSamSize)
library(tidyverse)
N <- seq(3*10^5,3*10^5,10000)
e <- 0.01
sample_01 <- N/(1+N*e^2)
#cor_sample <- correction(ns = sample_01,N = 3*10^5)
png(file="result_fig2.png",width=7,height=6,units = "in",res = 300)
plot(N,sample_01,ylim = c(0,10000),"o",xlim = c(10^3,10^6),col = "black",xlab = "population size N",
     ylab = "sample size" ,lwd = 1,lty = 5,pch = 1,cex = 0.2,yaxs="i")
lines(N,sample_02,col = "#e41a1c",lwd = 1,type = "o",lty = 5,pch = 2,cex = 0.2)
title("Simplified sample size estimation",adj = 0)
legend("topright", inset = c(0,0.03),
       legend = c("error margin = 0.01", "error margin = 0.02"),
       col = c('black', 
               "#e41a1c"),pch =c(1,2),
       pt.cex = 0.5, lwd = 0.8,
       cex = 1, 
       text.col = "black", 
       horiz = F )
dev.off()



#Evaluate for a grid of hypothesized proportion
correction <- function(ns,N){
  y = (ns * N)/(ns+N-1)
  return(y)
}
p.grid <- seq(0.5,0.5,0.01)
cissfuns <- list(ciss.wald, ciss.agresticoull,ciss.wilson)
ns <- sapply(p.grid, function(p) {
  unlist(lapply(cissfuns, function(f) f(p, d=0.01, alpha=0.01)))
})

ns_cor <- apply(ns,c(1,2),correction,N = 3*10^5)

png(file="result_fig1d.png",width=7,height=6,units = "in",res = 300)
matplot(p.grid, t(ns),type="l",xlab=expression(p[0]),ylab="n",lwd=2, xaxs = "i",yaxs = "i")
legend(x="topleft", c("Wald", "Wilson","Agresti-Coull"), col=1:3, lty=1:3,lwd=2)
title("CL = 99% and interval of \u00b1 0.01",adj = 0)
dev.off()
ns_df_01_01 <- as.data.frame(t(ns))



ns_cor_df_01_01 <- as.data.frame((apply(ns_df_01_01,c(1,2),correction,N = 10^5)))
matplot(p.grid, t(ns_cor),type="l",xlab=expression(p[0]),ylab="n",lwd=2)
legend(x="topleft", c("Wald", "Wilson","Agresti-Coull"), col=1:3, lty=1:3,lwd=2)


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


x <- wilson(p = 25/50, 0.95, 50)
x[[1]]
x[[1]]-x[[2]]
x[[1]]+x[[2]]
