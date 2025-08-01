


###计算F值和R2(无eaf版)
get_f_noeaf <- function(dat, F_value = 10) {
  if(is.null(dat$beta.exposure[1]) || is.na(dat$beta.exposure[1])) {
    print("数据不包含beta，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$se.exposure[1]) || is.na(dat$se.exposure[1])) {
    print("数据不包含se，无法计算F统计量")
    return(dat)
  }
  if(is.null(dat$samplesize.exposure[1]) || is.na(dat$samplesize.exposure[1])) {
    print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)
  }
  
  R2 <- (dat$beta.exposure^2) / ((dat$se.exposure^2 * dat$samplesize.exposure) + dat$beta.exposure^2)
  F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
  dat$R2 <- R2
  dat$F <- F
  dat <- subset(dat, F > F_value)
  return(dat)
}

#######计算F值和R2#######
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("数据不包含eaf，无法计算F统计量")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)}
  
  
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}




###steiger检验

steiger_test <- function(dat) {
  dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                    dat$se.exposure,
                                    dat$samplesize.exposure)
  dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                   dat$se.outcome,
                                   dat$samplesize.outcome)
  res_steiger <- mr_steiger(
    p_exp = dat$pval.exposure,
    p_out = dat$pval.outcome,
    n_exp = dat$samplesize.exposure,
    n_out = dat$samplesize.outcome,
    r_exp = dat$r.exposure,
    r_out = dat$r.outcome
  )
  res_steiger <- directionality_test(dat)
  
  return(res_steiger)
}


###power计算
results_binary <- function(N, alpha, R2xz, K, OR, epower) {
  threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
  f.value <- 1 + N * R2xz / (1 - R2xz)
  
  if (is.na(epower)) {
    
    b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
    
    v_MR <- (K * (1-K) - b_MR^2) / (N*R2xz)
    NCP <- b_MR^2 / v_MR
    
    # 2-sided test
    power <- 1 - pchisq(threschi, 1, NCP)
    data.frame(Parameter = c("Power", "NCP", "F-statistic"), Value = c(power, NCP, f.value), Description = c("", "Non-Centrality-Parameter", "The strength of the instrument"))    
    
    
  } else {
    
    # Calculation of sample size given power
    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(epower)
    Z  <- (z1 + z2)^2
    
    b_01 <- K * ( OR/ (1 + K * (OR - 1)) -1)
    f <- K * (1-K) - b_01^2
    N1 <- Z * f / (b_01^2 * R2xz)
    N1 <- ceiling(N1)
    data.frame(Parameter = "Sample Size", Value = N1)
    
  }
}


#####MR-PRESSO######
# PRESSO
mr_Presso<-function(dat,num=10000){
  library(TwoSampleMR)
  library(MRPRESSO)
  library(dplyr)
  set.seed(123)
  try (mr_presso_res<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                SignifThreshold = 0.05,NbDistribution = num))
  return(mr_presso_res)
  
}
mr_presso_pval<-function(mr_presso_res){ 
  try ( mr_presso_main<-mr_presso_res$`Main MR results`)
  try ( mr_presso_main[3:5,]<-NA) 
  return(mr_presso_main)
}


mr_presso_snp<-function(mr_presso_res,mr_presso_main,dat,type="list"){
  data_re<-list()
  if(type=="list"){
    for(i in 1:length(mr_presso_res)){
      res<-mr_presso_res[[i]]
      main<-mr_presso_main[[i]]
      data<-dat[[i]]
      try(if(is.na(main[2,6])==FALSE){
        outliers<-which(res$`Outlier Test`$Pvalue<0.05)
        data$mr_keep[outliers]<-FALSE
      })
      data_re[[i]]<-data
      names(data_re)[[i]]<-names(dat)[[i]]
    }
    return(data_re)
  }
  
  if(type=="data"){
    res<-mr_presso_res$`MR-PRESSO results`
    main<-mr_presso_main
    data<-dat
    try(if(is.na(main[2,6])==FALSE){
      outliers<-which(res$`Outlier Test`$Pvalue<0.05)
      data$mr_keep[outliers]<-FALSE
    })
    return(data)
  }
}



######自动选择IVW随即固定效应模型#####
choose_MR <- function(dat = dat) {
  res_hete <- NULL  
  if (nrow(dat) < 3) {
    res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
  } else {
    res_hete <- mr_heterogeneity(dat)
    if (res_hete$Q_pval[2] < 0.05) {
      res <- mr(dat, method_list = c(
        "mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_weighted_mode", "mr_simple_mode"
      ))
    } else {
      res <- mr(dat, method_list = c(
        "mr_egger_regression", "mr_weighted_median", "mr_ivw_fe", "mr_weighted_mode", "mr_simple_mode"
      ))
    }
  }
  AAA <- list(res_hete, res)
  return(list(AAA))
}



