eval_performance <- function(true_labels,predicted_labels,conf_int=0.95,method='Wald'){
  if (conf_int==0.9){
    z <- 1.645
  }else if (conf_int==0.95){
    z <- 1.96
  }else if (conf_int==0.98){
    z <- 2.33
  }else if (conf_int==0.99){
    z <- 2.576
  }else {
    stop("Confidence Intervals must be one of the real numbers: 0.90,0.95,0.98,0.99")
  }
  
  if (method=='Wald'){
    find_with_intervs <- function(Y,n,z){
      pi_hat <- Y/n
      pi_low <- pi_hat-z*sqrt(pi_hat*(1-pi_hat)/n)
      pi_upper <- pi_hat+z*sqrt(pi_hat*(1-pi_hat)/n)
      return(c(pi_low,pi_hat,pi_upper))
    }
  }else if (method=='Agresti-Coull'){
    find_with_intervs <- function(Y,n,z){
      pi_tilda <- (Y+0.5*z^2)/(n+z^2)
      n_tilda <- n+z^2
      pi_low <- pi_tilda-z*sqrt(pi_tilda*(1-pi_tilda)/n_tilda)
      pi_upper <- pi_tilda+z*sqrt(pi_tilda*(1-pi_tilda)/n_tilda)
      return(c(pi_low,pi_tilda,pi_upper))
    }
  }else if (method=='Wilson') {
    find_with_intervs <- function(Y,n,z){
      pi_hat <- Y/n
      pi_tilda <- (Y+0.5*z^2)/(n+z^2)
      n_tilda <- n+z^2
      pi_low <- pi_tilda-(z/n_tilda)*sqrt(pi_hat*(1-pi_hat)*n+0.25*z^2)
      pi_upper <- pi_tilda+(z/n_tilda)*sqrt(pi_hat*(1-pi_hat)*n+0.25*z^2)
      return(c(pi_low,pi_tilda,pi_upper))
    }
  }else{
    stop("Method must be one of the following: Wald,Agresti-Coull,Wilson")
  }
  df <- data.frame(true_labels = true_labels, predicted_labels = predicted_labels)
  
  tp <- df %>% filter(predicted_labels == "Positive") %>% filter(true_labels=="Positive")
  tp <- nrow(tp)
  
  fp <- df %>% filter(predicted_labels == "Positive") %>% filter(true_labels=="Negative")
  fp <- nrow(fp)
  
  tn <- df %>% filter(predicted_labels == "Negative") %>% filter(true_labels=="Negative")
  tn <- nrow(tn)
  
  fn <- df %>% filter(predicted_labels == "Negative") %>% filter(true_labels=="Positive")
  fn <- nrow(fn)
  
  sens <- find_with_intervs(Y=tp,n=(tp+fn),z=z)
  spec <- find_with_intervs(Y=tn,n=(tn+fp),z=z)
  ppv <- find_with_intervs(Y=tp,n=(tp+fp),z=z)
  npv <- find_with_intervs(Y=tn,n=(tn+fn),z=z)
  acc <- find_with_intervs(Y=(tn+tp),n=(tn+fn+tp+fp),z=z)
  #sens <- tp/(tp+fn)
  #spec <- tn/(tn+fp)
  
  #ppv <- tp/(tp+fp)
  #npv <- tn/(tn+fn)
  #acc <- (tn+tp)/(tn+fn+tp+fp)
  
  performance <- as.data.frame(cbind(sens[1],sens[2],sens[3],spec[1],spec[2],spec[3],
                       ppv[1],ppv[2],ppv[3],npv[1],npv[2],npv[3],acc[1],acc[2],acc[3],fp,fn))
  #performance <- t(performance)
  colnames(performance) <-  c('sensitivity_lower','sensitivity','sensitivity_upper','specificity_lower','specificity','specificity_upper',
                              'PPV_lower','PPV','PPV_upper','NPV_lower','NPV','NPV_upper','accuracy_lower','accuracy','accuracy_upper','fp','fn')
  #rownames(performance) <- c('sensitivity','specificity','PPV','NPV','accuracy')
  return(performance)
}
