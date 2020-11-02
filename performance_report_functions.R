calc_thresh <- function(df,outlier){
  #first is always Samples and second is always labels
  n <- ncol(df)
  npred <- n-2
  limits <- matrix(0,nrow = npred,ncol = 4)
  limits <- as.data.frame(limits)
  colnames(limits) <- c("predictor","lower","thresh","upper")
  for (i in 1:npred) {
    meas <- as.numeric(df[,i+2])
    meas_neg <- meas[df$labels=="Negative"]
    m <- mean(meas_neg)
    sd <- sd(meas_neg)
    thresh <- m + 3*sd
    id <- which(meas[df$labels=="Negative"]>(outlier)*thresh)
    if (length(id) != 0) {
      neg <- which(df$labels=="Negative")
      neg <- neg[-id]
      thresh <- mean(meas[neg]) + 3* sd(meas[neg])
      m <- mean(meas[neg])
      sd <- sd(meas[neg])
    }
    limits[i,2] <- m + 2*sd
    limits[i,3] <- m + 3*sd
    limits[i,4] <- m + 4*sd
  }
  limits$predictor <- as.character(colnames(df)[3:n])
  return(limits)
}

find_with_intervs <- function(Y,n,z){
  pi_hat <- Y/n
  pi_tilda <- (Y+0.5*z^2)/(n+z^2)
  n_tilda <- n+z^2
  pi_low <- pi_tilda-(z/n_tilda)*sqrt(pi_hat*(1-pi_hat)*n+0.25*z^2)
  pi_upper <- pi_tilda+(z/n_tilda)*sqrt(pi_hat*(1-pi_hat)*n+0.25*z^2)
  return(c(pi_low,pi_hat,pi_upper))}
calc_neg_thresh_wo_outliers <- function(df,iqr,thresh_sd){
  #first is always Samples and second is always labels
  n <- ncol(df)
  npred <- n-2
  limits <- matrix(0,nrow = npred,ncol = 7)
  limits <- as.data.frame(limits)
  colnames(limits) <- c("predictor","lower","thresh","upper","mean","sd","5sd")
  for (i in 1:npred) {
    meas <- as.numeric(df[,i+2])
    meas_neg <- meas[df$labels=="Negative"]
    box <- boxplot(meas_neg)
    lower_pos <- box$stats[2]
    upper_pos <- box$stats[4]
    iqr_pos <- iqr*(upper_pos - lower_pos)
    outlier_up <- upper_pos + iqr_pos
    #m <- mean(meas_neg[meas_neg<outlier_up])
    m <- mean(meas_neg)
    #sd <- sd(meas_neg[meas_neg<outlier_up])
    sd <- sd(meas_neg)
    thresh <- m + thresh_sd*sd
    limits[i,2] <- m + 2*sd
    limits[i,3] <- thresh
    limits[i,4] <- m + 4*sd
    limits[i,5] <- m
    limits[i,6] <- sd
    limits[i,7] <- m+5*sd
  }
  limits$predictor <- as.character(colnames(df)[3:n])
  return(limits)
}

calc_truncated_neg_thresh_wo_outliers <- function(df,iqr,thresh_sd){
  #first is always Samples and second is always labels
  n <- ncol(df)
  npred <- n-2
  limits <- matrix(0,nrow = npred,ncol = 7)
  limits <- as.data.frame(limits)
  colnames(limits) <- c("predictor","lower","thresh","upper","mean","sd","5sd")
  for (i in 1:npred) {
    meas <- as.numeric(df[,i+2])
    meas_neg <- meas[df$labels=="Negative"]
    box <- boxplot(meas_neg)
    lower_pos <- box$stats[2]
    upper_pos <- box$stats[4]
    iqr_pos <- iqr*(upper_pos - lower_pos)
    outlier_up <- upper_pos + iqr_pos
    m <- mean(meas_neg[meas_neg<outlier_up])
    #m <- mean(meas_neg)
    sd <- sd(meas_neg[meas_neg<outlier_up])
    #sd <- sd(meas_neg)
    truncated_meas <- meas_neg[meas_neg>=m]
    truncated_meas <- truncated_meas[truncated_meas<outlier_up]
    #h <- hist(truncated_meas)
    ##correct bins
    #if (h$counts[1]<0.85*max(h$counts)){
    #ind <- 2
    #if (h$counts[2]<0.85*max(h$counts)){
    #ind <- 3
    #}
    #tran <- h$breaks[ind]
    #truncated_meas <- truncated_meas[truncated_meas>tran]
    #}
    mirror_num <- m
    v <- truncated_meas-mirror_num
    symetric_v <- mirror_num-v
    modeled_gauss <- c(symetric_v,truncated_meas)
    #m <- mean(modeled_gauss)
    sd <- sd(modeled_gauss)
    #hist(modeled_gauss)
    thresh <- m + thresh_sd*sd
    limits[i,2] <- m + 2*sd
    limits[i,3] <- thresh
    limits[i,4] <- m + 4*sd
    limits[i,5] <- m
    limits[i,6] <- sd
    limits[i,7] <- m+5*sd
  }
  limits$predictor <- as.character(colnames(df)[3:n])
  return(limits)
}

sens_spec <- function(meas,labels,thresh){
  df <- data.frame(meas = meas, labels = labels)
  df$meas <- df$meas >= thresh
  df$meas <- df$meas + 0
  df <- df %>% mutate(meas = if_else(meas == 1,true = "Positive",false = "Negative"))
  
  tp <- df %>% filter(meas == "Positive") %>% filter(labels=="Positive")
  tp <- nrow(tp)
  
  fp <- df %>% filter(meas == "Positive") %>% filter(labels=="Negative")
  fp <- nrow(fp)
  
  tn <- df %>% filter(meas == "Negative") %>% filter(labels=="Negative")
  tn <- nrow(tn)
  
  fn <- df %>% filter(meas == "Negative") %>% filter(labels=="Positive")
  fn <- nrow(fn)
  
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  return(c(sens,spec,1-spec,fn,fp))
}


ciROCcustom <- function(rocit_emp, level){
  argClass <- class(rocit_emp)
  if(argClass != "rocit"){
    stop("Argument is not of class \"rocit\" ")
  }
  
  argMethod <- rocit_emp$method
  if(argMethod != "empirical"){
    stop("Rocit object method is not \"empirical\" ")
  }
  
  # initialize values
  pos_count <- rocit_emp$pos_count
  neg_count <- rocit_emp$neg_count
  pos_D <- rocit_emp$pos_D
  neg_D <- rocit_emp$neg_D
  TPR <- rocit_emp$TPR
  FPR <- rocit_emp$FPR
  c <- rocit_emp$Cutoff
  max_nD <- max(neg_D)
  min_nD <- min(neg_D)
  # approximate distributions of diagnostic in two groups
  ppdf_temp <- approxfun(density(pos_D))
  npdf_temp <- approxfun(density(neg_D))
  ppdf <- function(x){
    ifelse(is.na(ppdf_temp(x)), 0, ppdf_temp(x))
  }
  npdf <- function(x){
    ifelse(is.na(npdf_temp(x)), 0, npdf_temp(x))
  }
  
  
  # CDF for negative diagnostic
  ncdf_2 <- function(y){
    retval <- rep(NA, length(y))
    for(i in 1:length(y)){
      x <- y[i]
      if(x>max(neg_D)){retval[i] <- 1}else{
        if(x<min(neg_D)){retval[i] <- 0}else{
          bw <- (x-min(neg_D))/1000
          myX <- seq(min(neg_D), x, bw)
          myY <- npdf(myX)
          retval[i] <- sum(diff(myX)*(myY[-1]+myY[-length(myY)]))/2
        }
      }
    }
    return(retval)
  }
  # survival function for negative diagnostic
  nSurv <- function(x) {1-ncdf_2(x)}
  # groundwork for c_star
  binwidht <- (max_nD-min_nD)/10000
  DummyX <- seq(min_nD-0.001,max_nD+0.001,binwidht)
  DummyY <- nSurv(DummyX)
  c_starfun=approxfun(DummyY,DummyX)
  c_starfun2 <- function(x){
    ifelse(x==0, max_nD, ifelse(x==1,min_nD,c_starfun(x)))
  }
  # calculations
  c_star <- c_starfun2(FPR)
  #var_term2 <- (npdf(c_star)/ppdf(c_star))^2 * FPR * (1-FPR)/neg_count
  var_term1 <- TPR * (1-TPR)/pos_count
  SE_TPR <- sqrt(var_term1 )
  multiplier <- qnorm((1+level)/2)
  upper <- TPR + multiplier * SE_TPR
  lower <- TPR - multiplier * SE_TPR
  # cap CI limits within [0,1]
  upFun <- function(x) {min(c(1,x))}
  lowFun <- function(x) {max(c(0,x))}
  upper <- sapply(upper, upFun)
  lower <- sapply(lower, lowFun)
  # return
  TPR0 <- which(TPR == 0)
  TPR1 <- which(TPR == 1)
  lower[TPR1] <- 1
  upper[TPR0] <- 0
  returnval <- list(`ROC estimation method` = rocit_emp$method,
                    `Confidence level` = paste0(100*level,"%"),
                    FPR = FPR, TPR = TPR,
                    LowerTPR = lower, UpperTPR = upper,se = SE_TPR, mult = multiplier)
  return(returnval)
}

roc_anal <- function(ncov, n = 1000, output_dir,conc,type,w) {
  npred <- ncol(ncov)-2
  df <- matrix(0,nrow = npred,ncol = 10)
  df <- as.data.frame(df)
  colnames(df) <- c("predictor","lower","auc","upper","roc_thresh","sens","spec","fn","fp","youden")
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    minimum <- min(meas)
    maximum <- max(meas)
    steps <- seq(minimum,maximum,length.out = n)
    points <- matrix(10,nrow = length(steps),ncol = 5)
    for (j in 1:length(steps)) {
      points[j,] <- sens_spec(meas = meas,labels = ncov$labels,thresh = steps[j])
    }
    points <- as.data.frame(points)
    colnames(points) <- c("sens","spec","one_spec","fn","fp")
    points$id <- 1:nrow(points)
    points$thresh <- steps
    w_sens <- w[i-2]
    w_spec <- 1-w_sens
    points$y <- 2*(w_sens*points$sens+w_spec*points$spec)-1
    roc_thresh <- points$thresh[which(points$y==max(points$y))]
    if (length(roc_thresh)>1) {
      roc_thresh <- roc_thresh[order(roc_thresh,decreasing = F)]
      id_thresh <- round(length(roc_thresh)/2)
      roc_thresh <- roc_thresh[id_thresh]
    }
    write.csv(points,paste0(output_dir,"/",type,"_",colnames(ncov)[i],"_",conc,"_roc_data.csv"))
    library(ROCit)
    roc_empirical <- rocit(score = meas, class = as.factor(ncov$labels),
                           negref = "Negative")
    auc <- ciAUC(roc_empirical)
    roc_ci <- ciROCcustom(roc_empirical, level = 0.95)
    png(filename = paste0(output_dir,"/",type,"_",colnames(ncov)[i],"_",conc,".png"),width = 9,height = 6,units = "in",res = 300)
    plot(roc_empirical,values=F,YIndex = F,main = paste0(type,"-",colnames(ncov)[i],"-",conc))
    lines(roc_ci$FPR,roc_ci$LowerTPR,col= "red",lty =2)
    lines(roc_ci$FPR,roc_ci$UpperTPR,col= "red",lty =2)
    dev.off()
    df$lower[i-2] <- auc$lower
    df$auc[i-2] <- auc$AUC
    df$upper[i-2] <- auc$upper
    df$roc_thresh[i-2] <- roc_thresh
    df$predictor[i-2] <- colnames(ncov)[i]
    df$sens[i-2] <- points$sens[which(points$thresh == roc_thresh)]
    df$spec[i-2] <- points$spec[which(points$thresh == roc_thresh)]
    df$fn[i-2] <- points$fn[which(points$thresh == roc_thresh)]
    df$fp[i-2] <- points$fp[which(points$thresh == roc_thresh)]
    df$youden[i-2] <- points$y[which(points$thresh == roc_thresh)]
  }
  write.csv(df,paste0(output_dir,"/",type,"_thresh_performance.csv"))
  return(df)
}


plot_scatter_bar <- function(ncov,roc_perf,cutoff,output_dir,type,conc){
  gscatters <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    thresh1 <- cutoff$thresh[i-2]
    roc_thresh <- roc_perf$roc_thresh[i-2]
    general <- ggplot(ncov,aes(as.character(Sample),meas))+geom_point(aes(colour = labels),size = 0.5)+xlab("Samples")+ylab("Fluorescent Intensity") + 
      geom_hline(aes(yintercept=thresh1, linetype="cut-off"), color = "red")+
      geom_hline(aes(yintercept=roc_thresh, linetype="roc-threshold"), color = "blue")+
      scale_linetype_manual(name = "cutoff", values = c(1,2), 
                            guide = guide_legend(override.aes = list(color = c("red","blue"))))+
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc[i]))+
      theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0.2,size=2))
    
    #annotate(geom="text", x=10, y=(-5000), label="Negatives",color="black")+
    #annotate(geom="text", x=10, y=(thresh1+10000), label="Positives",color="black")+
    
    gen <- ggplot_build(general)
    general <- general + scale_y_continuous(breaks = sort(c(gen$layout$coord$labels(gen$layout$panel_params)[[1]]$y.major_source,round(thresh1,3),round(roc_thresh,3))))
    gscatters[[k]] <- ggplotGrob(general)
    general <- NULL
    gen <- NULL
    k <- k+1
    #print(general)
    #dev.off()
  }
  png(filename = paste0(output_dir,"/",type,"_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gscatters,nrow=1,top="Antigen Performance Analysis")
  dev.off()
  
  
  gbars <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    thresh1 <- cutoff$thresh[i-2]
    roc_thresh <- roc_perf$roc_thresh[i-2]
    df <- ncov[,c(1,i,2)]
    barplot <- ggplot(df,aes(x=factor(reorder(Sample,-df[,2])), y = df[,2],width=0.5,fill = labels)) +    
      geom_bar(stat="identity",position=position_nudge((x=0))) +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0,size=1))+
      xlab("Samples") + ylab("Fluorescent Intensity") + 
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc[i]))+ 
      geom_hline(aes(yintercept=thresh1, linetype="cut-off"), color = "red")+
      geom_hline(aes(yintercept=roc_thresh, linetype="roc-threshold"), color = "blue")+
      scale_linetype_manual(name = "cutoff", values = c(1,2), 
                            guide = guide_legend(override.aes = list(color = c("red","blue"))))
    bar <- ggplot_build(barplot)
    barplot <- barplot + scale_y_continuous(breaks = sort(c(bar$layout$coord$labels(bar$layout$panel_params)[[1]]$y.major_source,round(thresh1,3),round(roc_thresh,3))))
    #jpeg(file=paste0("igg_iga_igm_results/ig_all_bar_",colnames(igg)[i],".jpeg"),width=9,height=6,units = "in",quality=150,res=300)
    #print(barplot)
    #dev.off()
    gbars[[k]] <- ggplotGrob(barplot)
    k <- k+1
    barplot <- NULL
    bar <- NULL
  }
  
  png(filename = paste0(output_dir,"/",type,"_barplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gbars,nrow=1,top="Antigen Performance Analysis")
  dev.off()
}

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
      return(c(pi_low,pi_hat,pi_upper))
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


eval_performance_roc <- function(ncov, roc_perf, type, conc = "7.5ug",output_dir,normal = "raw"){
  npred <- ncol(ncov)-2
  predictions <- NULL
  predictions <- ncov[,c(1,2)]
  predictions <- as.data.frame(predictions)
  for (i in 1:npred) {
    meas <- as.numeric(ncov[,i+2])
    pred <- meas >= roc_perf$roc_thresh[i]
    pred <- pred + 0
    flag1 <- (abs(meas-roc_perf$roc_thresh[i])/roc_perf$roc_thresh[i])<0.2
    flag1 <- flag1 + 0
    predictions <- cbind(predictions,pred,flag1)
  }
  
  s <- seq(1,npred*2,2)
  #naming the predictions
  for (i in 1:length(s)) {
    
    colnames(predictions)[s[i]+2] <- paste0(colnames(ncov)[i+2],"_pred")
    colnames(predictions)[s[i]+3] <- paste0(colnames(ncov)[i+2],"_flag1")
    
  }
  # calc predictions
  predictions <- predictions %>% mutate(rule = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_pred == 1 )),
                                                       true = "Positive",false = "Negative")) %>%
    mutate(RBD = if_else(RBD_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1 = if_else(S1_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1_rbd_rule = if_else(condition = (RBD_pred == 1 & S1_pred == 1),
                                                       true = "Positive",false = "Negative"))
  
  predictions$rule_flag <- "empty"
  for (i in 1:nrow(predictions)) {
    if (predictions$RBD_flag1[i]==0) {
      if (predictions$RBD_pred[i]==1){
        predictions$rule_flag[i] <- "Positive"
      }
      if (predictions$RBD_pred[i]==0){
        predictions$rule_flag[i] <- "Negative"
      }
    }
    if (predictions$RBD_flag1[i]==1) {
      if (predictions$RBD_pred[i]==1 & (predictions$N_pred[i] == 1 | predictions$S1_pred[i] == 1)){
        predictions$rule_flag[i] <- "Positive"
      } else {
        predictions$rule_flag[i] <- "Negative"
      }
    }
  }
  
  performance_rule <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule,conf_int = 0.95,method = "Wald")
  performance_rule$rule <- "rule"
  performance_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$RBD,conf_int = 0.95,method = "Wald")
  performance_rbd$rule <- "RBD only"
  performance_rule_flag <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_flag,conf_int = 0.95,method = "Wald")
  performance_rule_flag$rule <- "rule flag"
  performance_s1 <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1,conf_int = 0.95,method = "Wald")
  performance_s1$rule <- "S1 only"
  
  performance_s1_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1_rbd_rule,conf_int = 0.95,method = "Wald")
  performance_s1_rbd$rule <- "RBD + S1"
  
  
  performance <- bind_rows(performance_rule,performance_rbd,performance_rule_flag,performance_s1,performance_s1_rbd)
  performance$conc <- conc
  performance$thresh_type <- "roc_thresh"
  performance$thresh_rbd <- roc_perf$roc_thresh[roc_perf$predictor=="RBD"]
  performance$thresh_s1 <- roc_perf$roc_thresh[roc_perf$predictor=="S1"]
  performance$thresh_n <- roc_perf$roc_thresh[roc_perf$predictor=="N"]
  performance$auc_RBD_lower <- roc_perf$lower[roc_perf$predictor=="RBD"]
  performance$auc_RBD <- roc_perf$auc[roc_perf$predictor=="RBD"]
  performance$auc_RBD_upper <- roc_perf$upper[roc_perf$predictor=="RBD"]
  performance$type <- type
  performance$normalization <- normal
  
  write.csv(performance,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_performance_df.csv"))
  write.csv(predictions,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_predictions_df.csv"))
  return(performance)
}

eval_performance_thresh <- function(ncov, cutoff, roc_perf, type, conc = "7.5ug",output_dir,normal = "raw"){
  npred <- ncol(ncov)-2
  predictions <- NULL
  predictions <- ncov[,c(1,2)]
  predictions <- as.data.frame(predictions)
  for (i in 1:npred) {
    meas <- as.numeric(ncov[,i+2])
    pred <- meas >= cutoff$thresh[i]
    pred <- pred + 0
    flag1 <- (meas>=cutoff$lower[i] & meas <= cutoff$upper[i])
    flag1 <- flag1 + 0
    predictions <- cbind(predictions,pred,flag1)
  }
  
  s <- seq(1,npred*2,2)
  #naming the predictions
  for (i in 1:length(s)) {
    
    colnames(predictions)[s[i]+2] <- paste0(colnames(ncov)[i+2],"_pred")
    colnames(predictions)[s[i]+3] <- paste0(colnames(ncov)[i+2],"_flag1")
    
  }
  # calc predictions
  predictions <- predictions %>% mutate(rule = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_pred == 1 )),
                                                       true = "Positive",false = "Negative")) %>%
    mutate(RBD = if_else(RBD_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1 = if_else(S1_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1_rbd_rule = if_else(condition = (RBD_pred == 1 & S1_pred == 1),
                                                              true = "Positive",false = "Negative"))
  
  predictions$rule_flag <- "empty"
  for (i in 1:nrow(predictions)) {
    if (predictions$RBD_flag1[i]==0) {
      if (predictions$RBD_pred[i]==1){
        predictions$rule_flag[i] <- "Positive"
      }
      if (predictions$RBD_pred[i]==0){
        predictions$rule_flag[i] <- "Negative"
      }
    }
    if (predictions$RBD_flag1[i]==1) {
      if (predictions$RBD_pred[i]==1 & (predictions$N_pred[i] == 1 | predictions$S1_pred[i] == 1)){
        predictions$rule_flag[i] <- "Positive"
      } else {
        predictions$rule_flag[i] <- "Negative"
      }
    }
  }
  
  performance_rule <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule,conf_int = 0.95,method = "Wald")
  performance_rule$rule <- "rule"
  performance_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$RBD,conf_int = 0.95,method = "Wald")
  performance_rbd$rule <- "RBD only"
  performance_rule_flag <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_flag,conf_int = 0.95,method = "Wald")
  performance_rule_flag$rule <- "rule flag"
  
  performance_s1 <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1,conf_int = 0.95,method = "Wald")
  performance_s1$rule <- "S1 only"
  
  performance_s1_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1_rbd_rule,conf_int = 0.95,method = "Wald")
  performance_s1_rbd$rule <- "RBD + S1"
  
  performance <- bind_rows(performance_rule,performance_rbd,performance_rule_flag,performance_s1,performance_s1_rbd)
  
  performance$conc <- conc
  performance$thresh_type <- "negative_thresh"
  performance$thresh_rbd <- cutoff$thresh[cutoff$predictor=="RBD"]
  performance$thresh_s1 <- cutoff$thresh[cutoff$predictor=="S1"]
  performance$thresh_n <- cutoff$thresh[cutoff$predictor=="N"]
  performance$auc_RBD_lower <- roc_perf$lower[roc_perf$predictor=="RBD"]
  performance$auc_RBD <- roc_perf$auc[roc_perf$predictor=="RBD"]
  performance$auc_RBD_upper <- roc_perf$upper[roc_perf$predictor=="RBD"]
  performance$type <- type
  performance$normalization <- normal
  
  write.csv(performance,paste0(output_dir,"/",type,"_",conc,"_neg_thresh_performance_df.csv"))
  write.csv(predictions,paste0(output_dir,"/",type,"_",conc,"_neg_thresh_predictions_df.csv"))
  return(performance)
}

eval_performance_roc_anal <- function(ncov, roc_perf, type, conc = "7.5ug",output_dir,normal = "raw",outlier_thresh){
  npred <- ncol(ncov)-2
  predictions <- NULL
  predictions <- ncov[,c(1,2)]
  predictions <- as.data.frame(predictions)
  for (i in 1:npred) {
    meas <- as.numeric(ncov[,i+2])
    pred <- meas >= roc_perf$roc_thresh[i]
    pred <- pred + 0
    flag1 <- (abs(meas-roc_perf$roc_thresh[i])/roc_perf$roc_thresh[i])<=outlier_thresh
    flag1 <- flag1 + 0
    predictions <- cbind(predictions,pred,flag1)
  }
  
  s <- seq(1,npred*2,2)
  #naming the predictions
  for (i in 1:length(s)) {
    
    colnames(predictions)[s[i]+2] <- paste0(colnames(ncov)[i+2],"_pred")
    colnames(predictions)[s[i]+3] <- paste0(colnames(ncov)[i+2],"_flag1")
    
  }
  # calc predictions
  predictions <- predictions %>% mutate(rule = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_pred == 1 )),
                                                       true = "Positive",false = "Negative")) %>%
    mutate(RBD = if_else(RBD_pred==1,true = "Positive",false = "Negative")) %>%
    mutate(N = if_else(N_pred==1,true = "Positive",false = "Negative")) %>%
    mutate(n_and_s = if_else(condition = (N_pred == 1 & S1_pred == 1),
                             true = "Positive",false = "Negative")) %>%
    mutate(rule_n = if_else(condition = (N_pred == 1 & (RBD_pred==1 | S1_pred == 1 )),
                          true = "Positive",false = "Negative")) %>%
    mutate(rule_s = if_else(condition = (S1_pred == 1 & (RBD_pred==1 | N_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rule_rbd_or = if_else(condition = (RBD_pred == 1 | (S1_pred==1 & N_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rule_s_or = if_else(condition = (S1_pred == 1 | (RBD_pred==1 & N_pred == 1 )),
                                 true = "Positive",false = "Negative")) %>%
    mutate(rule_n_or = if_else(condition = (N_pred == 1 | (RBD_pred==1 & S1_pred == 1 )),
                               true = "Positive",false = "Negative")) 
    
  
  predictions <- predictions %>% mutate(s1 = if_else(S1_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1_rbd_rule = if_else(condition = (RBD_pred == 1 & S1_pred == 1),
                                                              true = "Positive",false = "Negative"),
                                        n_rbd_rule = if_else(condition = (RBD_pred == 1 & N_pred == 1),
                                                             true = "Positive",false = "Negative"),
                                        all = if_else(condition = (RBD_pred == 1 & N_pred == 1 & S1_pred == 1),
                                                      true = "Positive",false = "Negative")) 
  
  predictions$rule_flag <- "empty"
  for (i in 1:nrow(predictions)) {
    if (predictions$RBD_flag1[i]==0) {
      if (predictions$RBD_pred[i]==1){
        predictions$rule_flag[i] <- "Positive"
      }
      if (predictions$RBD_pred[i]==0){
        predictions$rule_flag[i] <- "Negative"
      }
    }
    if (predictions$RBD_flag1[i]==1) {
      if (predictions$RBD_pred[i]==1 & (predictions$N_pred[i] == 1 | predictions$S1_pred[i] == 1)){
        predictions$rule_flag[i] <- "Positive"
      } else {
        predictions$rule_flag[i] <- "Negative"
      }
    }
  }
  
  performance_rule <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule,conf_int = 0.95,method = "Wald")
  performance_rule$rule <- "rule"
  
  performance_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$RBD,conf_int = 0.95,method = "Wald")
  performance_rbd$rule <- "RBD only"
  performance_N <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$N,conf_int = 0.95,method = "Wald")
  performance_N$rule <- "N only"
  performance_rule_flag <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_flag,conf_int = 0.95,method = "Wald")
  performance_rule_flag$rule <- "rule flag"
  performance_s1 <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1,conf_int = 0.95,method = "Wald")
  performance_s1$rule <- "S1 only"
  
  performance_all <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$all,conf_int = 0.95,method = "Wald")
  performance_all$rule <- "all"
  
  performance_s1_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1_rbd_rule,conf_int = 0.95,method = "Wald")
  performance_s1_rbd$rule <- "RBD + S1"
  
  performance_n_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$n_rbd_rule,conf_int = 0.95,method = "Wald")
  performance_n_rbd$rule <- "RBD + N"
  
  performance_n_and_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$n_and_s,conf_int = 0.95,method = "Wald")
  performance_n_and_s$rule <- "N + S"
  
  performance_rule_n <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_n,conf_int = 0.95,method = "Wald")
  performance_rule_n$rule <- "rule N"
  
  performance_rule_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_s,conf_int = 0.95,method = "Wald")
  performance_rule_s$rule <- "rule S"
  
  performance_rule_rbd_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_rbd_or,conf_int = 0.95,method = "Wald")
  performance_rule_rbd_or$rule <- "rule RBD OR"
  
  performance_rule_s_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_s_or,conf_int = 0.95,method = "Wald")
  performance_rule_s_or$rule <- "rule S OR"
  
  performance_rule_n_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_n_or,conf_int = 0.95,method = "Wald")
  performance_rule_n_or$rule <- "rule N OR"
  
  performance <- bind_rows(performance_rule,performance_rbd,performance_rule_flag,performance_s1,
                           performance_N,performance_s1_rbd, performance_n_rbd, performance_all,
                           performance_n_and_s, performance_rule_n, performance_rule_s, performance_rule_rbd_or,
                           performance_rule_s_or,performance_rule_n_or)
  performance$conc <- conc
  performance$thresh_type <- "custom"
  performance$thresh_rbd <- roc_perf$roc_thresh[roc_perf$predictor=="RBD"]
  performance$thresh_s1 <- roc_perf$roc_thresh[roc_perf$predictor=="S1"]
  performance$thresh_n <- roc_perf$roc_thresh[roc_perf$predictor=="N"]
  #performance$auc_RBD_lower <- roc_perf$lower[roc_perf$predictor=="RBD"]
  #performance$auc_RBD <- roc_perf$auc[roc_perf$predictor=="RBD"]
  #performance$auc_RBD_upper <- roc_perf$upper[roc_perf$predictor=="RBD"]
  performance$type <- type
  performance$normalization <- normal
  
  #write.csv(performance,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_performance_df.csv"))
  #write.csv(predictions,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_predictions_df.csv"))
  return(list(performance,predictions))
}

plot_scatter_bar_no_thresh <- function(ncov,output_dir,type,conc){
  gscatters <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    general <- ggplot(ncov,aes(as.character(Sample),meas))+geom_point(aes(colour = labels),size = 0.5)+xlab("Samples")+ylab("Fluorescent Intensity") + 
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc))+
      theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0.2,size=2))
    
    #annotate(geom="text", x=10, y=(-5000), label="Negatives",color="black")+
    #annotate(geom="text", x=10, y=(thresh1+10000), label="Positives",color="black")+
    
    gen <- ggplot_build(general)
    general <- general 
    gscatters[[k]] <- ggplotGrob(general)
    general <- NULL
    gen <- NULL
    k <- k+1
    #print(general)
    #dev.off()
  }
  png(filename = paste0(output_dir,"/",type,"_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gscatters,nrow=1,top="Antigen Performance Analysis")
  dev.off()
  
  
  gbars <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    df <- ncov[,c(1,i,2)]
    barplot <- ggplot(df,aes(x=factor(reorder(Sample,-df[,2])), y = df[,2],width=0.5,fill = labels)) +    
      geom_bar(stat="identity",position=position_nudge((x=0))) +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0,size=1))+
      xlab("Samples") + ylab("Fluorescent Intensity") + 
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc))
    bar <- ggplot_build(barplot)
    barplot <- barplot 
    #jpeg(file=paste0("igg_iga_igm_results/ig_all_bar_",colnames(igg)[i],".jpeg"),width=9,height=6,units = "in",quality=150,res=300)
    #print(barplot)
    #dev.off()
    gbars[[k]] <- ggplotGrob(barplot)
    k <- k+1
    barplot <- NULL
    bar <- NULL
  }
  
  png(filename = paste0(output_dir,"/",type,"_barplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gbars,nrow=1,top="Antigen Performance Analysis")
  dev.off()
}

eval_performance_neg_thresh <- function(ncov, thresh, type, conc = "7.5ug",output_dir,normal = "raw",outlier_thresh){
  npred <- ncol(ncov)-2
  predictions <- NULL
  predictions <- ncov[,c(1,2)]
  predictions <- as.data.frame(predictions)
  for (i in 1:npred) {
    meas <- as.numeric(ncov[,i+2])
    pred <- meas >= thresh$thresh[i]
    pred <- pred + 0
    flag1 <- (abs(meas-thresh$thresh[i])/thresh$thresh[i])<=outlier_thresh
    flag1 <- flag1 + 0
    predictions <- cbind(predictions,pred,flag1)
  }
  
  s <- seq(1,npred*2,2)
  #naming the predictions
  for (i in 1:length(s)) {
    
    colnames(predictions)[s[i]+2] <- paste0(colnames(ncov)[i+2],"_pred")
    colnames(predictions)[s[i]+3] <- paste0(colnames(ncov)[i+2],"_flag1")
    
  }
  # calc predictions
  predictions <- predictions %>% mutate(rule = if_else(condition = (RBD_pred == 1 & (N_pred==1 | S1_pred == 1 )),
                                                       true = "Positive",false = "Negative")) %>%
    mutate(RBD = if_else(RBD_pred==1,true = "Positive",false = "Negative")) %>%
    mutate(N = if_else(N_pred==1,true = "Positive",false = "Negative")) %>%
    mutate(n_and_s = if_else(condition = (N_pred == 1 & S1_pred == 1),
                             true = "Positive",false = "Negative")) %>%
    mutate(rule_n = if_else(condition = (N_pred == 1 & (RBD_pred==1 | S1_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rule_s = if_else(condition = (S1_pred == 1 & (RBD_pred==1 | N_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rule_rbd_or = if_else(condition = (RBD_pred == 1 | (S1_pred==1 & N_pred == 1 )),
                                 true = "Positive",false = "Negative")) %>%
    mutate(rule_s_or = if_else(condition = (S1_pred == 1 | (RBD_pred==1 & N_pred == 1 )),
                               true = "Positive",false = "Negative")) %>%
    mutate(rule_n_or = if_else(condition = (N_pred == 1 | (RBD_pred==1 & S1_pred == 1 )),
                               true = "Positive",false = "Negative")) %>%
    mutate(or_all = if_else(condition = (N_pred == 1 | (RBD_pred==1 | S1_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rbd_or_s = if_else(condition = ((RBD_pred==1 | S1_pred == 1 )),
                            true = "Positive",false = "Negative")) %>%
    mutate(rbd_or_n = if_else(condition = (N_pred == 1 | RBD_pred==1),
                            true = "Positive",false = "Negative")) %>%
    mutate(n_or_s = if_else(condition = (N_pred == 1 | S1_pred == 1 ),
                            true = "Positive",false = "Negative")) %>%
    mutate(maj = if_else(condition = (RBD_pred + N_pred + S1_pred)>=2,
                            true = "Positive",false = "Negative"))
  
  
  predictions <- predictions %>% mutate(s1 = if_else(S1_pred==1,true = "Positive",false = "Negative"))
  
  predictions <- predictions %>% mutate(s1_rbd_rule = if_else(condition = (RBD_pred == 1 & S1_pred == 1),
                                                              true = "Positive",false = "Negative"),
                                        n_rbd_rule = if_else(condition = (RBD_pred == 1 & N_pred == 1),
                                                             true = "Positive",false = "Negative"),
                                        all = if_else(condition = (RBD_pred & N_pred & S1_pred),
                                                      true = "Positive",false = "Negative")) 
  
  predictions$rule_flag <- "empty"
  for (i in 1:nrow(predictions)) {
    if (predictions$RBD_flag1[i]==0) {
      if (predictions$RBD_pred[i]==1){
        predictions$rule_flag[i] <- "Positive"
      }
      if (predictions$RBD_pred[i]==0){
        predictions$rule_flag[i] <- "Negative"
      }
    }
    if (predictions$RBD_flag1[i]==1) {
      if (predictions$RBD_pred[i]==1 & (predictions$N_pred[i] == 1 | predictions$S1_pred[i] == 1)){
        predictions$rule_flag[i] <- "Positive"
      } else {
        predictions$rule_flag[i] <- "Negative"
      }
    }
  }
  
  predictions$marginal_rule <- "Negative"
  
  for (i in 1:nrow(predictions)) {
    if (predictions$RBD_pred[i] == 1 & predictions$RBD_flag1[i]==0){
      if (predictions$N_pred[i] == 1 | predictions$S1_pred[i] == 1){
        predictions$marginal_rule[i] <- "Positive"
      }
      else if ((predictions$N_pred[i] == 0 & predictions$N_flag1[i] == 1) | (predictions$S1_pred[i] == 0 & predictions$S1_flag1[i] == 1)){
        predictions$marginal_rule[i] <- "Marginally Positive"
      }
    } else {
      if (predictions$N_pred[i] == 1 & predictions$N_flag1[i] == 0){
        if (predictions$RBD_pred[i] == 1 | predictions$S1_pred[i] == 1){
          predictions$marginal_rule[i] <- "Positive"
        } else {
          if ((predictions$RBD_pred[i] == 0 & predictions$RBD_flag1[i] == 1) | (predictions$S1_pred[i] == 0 & predictions$S1_flag1[i] == 1)){
            predictions$marginal_rule[i] <- "Marginally Positive"
          }
        }
      } else {
        if ((predictions$RBD_pred[i] == 1 & predictions$RBD_flag1[i]==1)&(predictions$N_pred[i] == 1 & predictions$N_flag1[i] == 1)&(predictions$S1_pred[i] == 1)){
          predictions$marginal_rule[i] <- "Positive"
        } else if ((predictions$RBD_pred[i] == 1 & predictions$N_pred[i] == 1)|(predictions$RBD_pred[i] == 1 & predictions$S1_pred[i] == 1)|(predictions$S1_pred[i] == 1 & predictions$N_pred[i] == 1)){
          predictions$marginal_rule[i] <- "Marginally Positive"
        }
      }
    }
  }
  
  predictions$marginal_positive <- "Negative"
  predictions$marginal_positive[grepl(pattern = "Positive",x = predictions$marginal_rule)] <- "Positive"
  
  predictions$marginal_negative <- predictions$marginal_rule
  predictions$marginal_negative[grepl(pattern = "Marginally",x = predictions$marginal_rule)] <- "Negative"
  
  performance_rule <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule,conf_int = 0.95,method = "Wilson")
  performance_rule$rule <- "rule"
  
  performance_marginal_pos <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$marginal_positive,conf_int = 0.95,method = "Wilson")
  performance_marginal_pos$rule <- "marginal positive"
  performance_marginal_neg <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$marginal_negative,conf_int = 0.95,method = "Wilson")
  performance_marginal_neg$rule <- "marginal negative"
  
  
  
  performance_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$RBD,conf_int = 0.95,method = "Wilson")
  performance_rbd$rule <- "RBD only"
  performance_N <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$N,conf_int = 0.95,method = "Wilson")
  performance_N$rule <- "N only"
  performance_rule_flag <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_flag,conf_int = 0.95,method = "Wilson")
  performance_rule_flag$rule <- "rule flag"
  performance_s1 <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1,conf_int = 0.95,method = "Wilson")
  performance_s1$rule <- "S1 only"
  
  performance_all <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$all,conf_int = 0.95,method = "Wilson")
  performance_all$rule <- "all"
  
  performance_s1_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$s1_rbd_rule,conf_int = 0.95,method = "Wilson")
  performance_s1_rbd$rule <- "RBD + S1"
  
  performance_n_rbd <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$n_rbd_rule,conf_int = 0.95,method = "Wilson")
  performance_n_rbd$rule <- "RBD + N"
  
  performance_n_and_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$n_and_s,conf_int = 0.95,method = "Wilson")
  performance_n_and_s$rule <- "N + S"
  
  performance_rule_n <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_n,conf_int = 0.95,method = "Wilson")
  performance_rule_n$rule <- "rule N"
  
  performance_rule_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_s,conf_int = 0.95,method = "Wilson")
  performance_rule_s$rule <- "rule S"
  
  performance_rule_rbd_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_rbd_or,conf_int = 0.95,method = "Wilson")
  performance_rule_rbd_or$rule <- "rule RBD OR"
  
  performance_rule_s_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_s_or,conf_int = 0.95,method = "Wilson")
  performance_rule_s_or$rule <- "rule S OR"
  
  performance_rule_n_or <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rule_n_or,conf_int = 0.95,method = "Wilson")
  performance_rule_n_or$rule <- "rule N OR"
  
  performance_n_or_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$n_or_s,conf_int = 0.95,method = "Wilson")
  performance_n_or_s$rule <- "N or S1"
  
  performance_rbd_or_s <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rbd_or_s,conf_int = 0.95,method = "Wilson")
  performance_rbd_or_s$rule <- "RBD or S1"
  
  performance_rbd_or_n <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$rbd_or_n,conf_int = 0.95,method = "Wilson")
  performance_rbd_or_n$rule <- "RBD or N"
  
  performance_or_all <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$or_all,conf_int = 0.95,method = "Wilson")
  performance_or_all$rule <- "ALL OR"
  
  performance_maj <- eval_performance(true_labels = ncov$labels,predicted_labels = predictions$maj,conf_int = 0.95,method = "Wilson")
  performance_maj$rule <- "Majority"
  
  performance <- bind_rows(performance_rule,performance_rbd,performance_rule_flag,performance_s1,
                           performance_N,performance_s1_rbd, performance_n_rbd, performance_all,
                           performance_n_and_s, performance_rule_n, performance_rule_s, performance_rule_rbd_or,
                           performance_rule_s_or,performance_rule_n_or,performance_n_or_s,performance_rbd_or_s,
                           performance_rbd_or_n,performance_or_all,performance_maj,performance_marginal_pos,
                           performance_marginal_neg)
  performance$conc <- conc
  performance$thresh_type <- "custom"
  performance$thresh_rbd <- thresh$thresh[thresh$predictor=="RBD"]
  performance$thresh_s1 <- thresh$thresh[thresh$predictor=="S1"]
  performance$thresh_n <- thresh$thresh[thresh$predictor=="N"]
  #performance$auc_RBD_lower <- roc_perf$lower[roc_perf$predictor=="RBD"]
  #performance$auc_RBD <- roc_perf$auc[roc_perf$predictor=="RBD"]
  #performance$auc_RBD_upper <- roc_perf$upper[roc_perf$predictor=="RBD"]
  performance$type <- type
  performance$normalization <- normal
  
  #write.csv(performance,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_performance_df.csv"))
  #write.csv(predictions,paste0(output_dir,"/",type,"_",conc,"_roc_thresh_predictions_df.csv"))
  return(list(performance,predictions))
}

plot_scatter_bar_neg_thresh <- function(ncov,cutoff,output_dir,type,conc){
  gscatters <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    thresh1 <- cutoff$thresh[i-2]
    general <- ggplot(ncov,aes(as.character(Sample),meas))+geom_point(aes(colour = labels),size = 0.5)+xlab("Samples")+ylab("Fluorescent Intensity") + 
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc))+
      geom_hline(yintercept=thresh1, linetype = "dashed",colour = "red")+
      theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0.2,size=2))
    
    #annotate(geom="text", x=10, y=(-5000), label="Negatives",color="black")+
    #annotate(geom="text", x=10, y=(thresh1+10000), label="Positives",color="black")+
    
    gen <- ggplot_build(general)
    general <- general 
    gscatters[[k]] <- ggplotGrob(general)
    general <- NULL
    gen <- NULL
    k <- k+1
    #print(general)
    #dev.off()
  }
  png(filename = paste0(output_dir,"/",type,"_scatterplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gscatters,nrow=1,top="Antigen Performance Analysis")
  dev.off()
  
  
  gbars <- NULL
  k <- 1
  for (i in 3:ncol(ncov)) {
    meas <- as.numeric(ncov[,i])
    thresh1 <- cutoff$thresh[i-2]
    df <- ncov[,c(1,i,2)]
    barplot <- ggplot(df,aes(x=factor(reorder(Sample,-df[,2])), y = df[,2],width=0.5,fill = labels)) +    
      geom_bar(stat="identity",position=position_nudge((x=0))) +theme(plot.title = element_text(hjust=0)) + theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust=0,size=1))+
      geom_hline(yintercept=thresh1, linetype = "dashed",colour = "red")+
      xlab("Samples") + ylab("Fluorescent Intensity") + 
      ggtitle(paste0(type,"-",colnames(ncov)[i],"-",conc))
    bar <- ggplot_build(barplot)
    barplot <- barplot 
    #jpeg(file=paste0("igg_iga_igm_results/ig_all_bar_",colnames(igg)[i],".jpeg"),width=9,height=6,units = "in",quality=150,res=300)
    #print(barplot)
    #dev.off()
    gbars[[k]] <- ggplotGrob(barplot)
    k <- k+1
    barplot <- NULL
    bar <- NULL
  }
  
  png(filename = paste0(output_dir,"/",type,"_barplot.png"),width = 18,height = 6,units = "in",res = 600)
  gridExtra::grid.arrange(grobs=gbars,nrow=1,top="Antigen Performance Analysis")
  dev.off()
}
