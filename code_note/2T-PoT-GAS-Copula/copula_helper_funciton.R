

calculate_marginal_label <- function( model, threshold){
  if( length(threshold) == 2 ){
    ret <- model$df.fitted$return
    
    threshold_0 <- threshold[["threshold_0"]]
    threshold_1 <- threshold[["threshold_1"]]
    #计算阈值
    # threshold_0 <- sort(ret[ret>0])[floor(length(ret[ret>0])*threshold_0)]
    # threshold_1 <- abs(sort(ret[ret<0])[floor(length(ret[ret<0])*(1-threshold_1))])
    threshold_0 <- sort(ret)[floor(length(ret)*threshold_0)]
    threshold_1 <- sort(ret)[floor(length(ret)*(1-threshold_1))]
    threshold <- c( threshold_0 = threshold_0, threshold_1  = threshold_1) 
    # 计算Pot_y
    ret = ret[ which.max( !is.na(model$tau) ):length( !is.na(model$tau) ) ]
    Pot_y = ret
    Pot_y[ret>0] = Pot_y[ret>0] - threshold_0
    Pot_y[ret<0] = Pot_y[ret<0] - threshold_1
    Pot_y[ret>0] = ifelse( Pot_y[ret>0]>0, Pot_y[ret>0], 0)
    Pot_y[ret<0] = ifelse( Pot_y[ret<0]<0, Pot_y[ret<0], 0)
  }else if( length(threshold) == 1 ){
    ret <- model$df.fitted$return
    #计算阈值
    
    threshold <- sort(ret)[floor(length(ret)*threshold)]
    # 计算Pot_y
    ret = ret[ which.max( !is.na(model$tau) ):length( !is.na(model$tau) ) ]
    Pot_y <- (ret-threshold)
    Pot_y = ifelse(Pot_y>0,Pot_y,0)
  }
 
  Pot_y_label <- ifelse( Pot_y == 0, "0", "1")
  return( list( Pot_y = Pot_y, label = Pot_y_label) )
}

calculate_Observation_driven <- function( Z, parameters, marginal_label){
  z1 <- Z[[1]]
  z2 <- Z[[2]]
  y <- parameters[1]
  a <- 1 + exp(y)
  # 定义a对似然函数求偏导
  l00 <- ( 
    -(z1^a+z2^a)^(1/a)
    *(-( (1/a)^2 )*log( z1^a+z2^a ) + (1/a)*( z1^a*log(z1) + z2^a*log(z2) )/(z1^a+z2^a) ) 
         ) 
  
  l10 <- (
    l00 - ( (1/a)^2 )*log( z1^a+z2^a ) + ( (1/a)-1 )*( z1^a*log(z1) + z2^a*log(z2) )/(z1^a+z2^a) + log(z1)
         )
  
  l01 <- (
    l00 - ( (1/a)^2 )*log( z1^a+z2^a ) + ( (1/a)-1 )*( z1^a*log(z1) + z2^a*log(z2) )/(z1^a+z2^a) + log(z2)
         )
  
  l11 <- ( 
    l00 + ( log(z1) + log(z2) ) - ( (1/a)^2 )*log( z1^a+z2^a )
           + ( (1/a)-2 )*( z1^a*log(z1) + z2^a*log(z2) )/(z1^a+z2^a) 
           + ( 1 - l00)/( (z1^a+z2^a)^(1/a) + a - 1 )
         )
  # if(is.na( get( paste0('l',marginal_label) ) )){
    # print( get( paste0('l',marginal_label) ) )
    # print( y )
    # print( a )
    # print( marginal_label )
    # print("数据包含NaN，程序终止。")
    # browser()
    # stop("数据包含NaN，程序终止。")
  # }
  return( get( paste0('l',marginal_label) )*(a-1) )
}

calculate_yt <- function( Z, parameters, marginal_label){
  beta0 <- parameters[1]
  beta1 <- parameters[2]
  beta2 <- parameters[3]
  
  y <- rep(0,length(Z[[1]]))
  u <- rep(0,length(Z[[1]]))
  y[1] <- beta0
  u[1] <- calculate_Observation_driven( c( Z[[1]][1], Z[[2]][1]), y[1], marginal_label[1])
  for (i in 2:length(Z[[1]])) {
    y[i] <- beta0 + beta1*y[i-1]+ beta2*u[i-1]
    u[i] <- calculate_Observation_driven( c( Z[[1]][i], Z[[2]][i]), y[i], marginal_label[i])
  }

  return(list(y,u))
}

llh_tvc <- function( Z, parameters, marginal_label){
  z1 <- Z[[1]]
  z2 <- Z[[2]]
  yt <- calculate_yt( Z, parameters, marginal_label)[[1]]
  at <- 1 + exp(yt)
  # # 定义基础似然函数
  # f00 <- ( 
  #   exp( -(z1^at+z2^at)^(1/at) )
  # ) 
  # 
  # f10 <- (
  #   f00*( z1^at + z2^at )^( 1/at-1 )*z1^( at-1 )
  # )
  # 
  # f01 <- (
  #   f00*( z1^at + z2^at )^( 1/at-1 )*z2^( at-1 )
  # )
  # 
  # f11 <- ( 
  #   f00*( ( z1^at + z2^at )^( 2/at-2 )*z1^( at-1 )*z2^( at-1 ) + ( z1^at + z2^at )^( 1/at-2 )*z1^( at-1 )*z2^( at-1 )*( at-1 ) )
  # )
  # 定义基础似然函数
  f00 <- ( 
    exp( -(z1^at+z2^at)^(1/at) )
  ) 
  
  f10 <- (
    f00*( z1^at + z2^at )^( 1/at-1 )*z1^( at-1 )*exp(z1)
  )
  
  f01 <- (
    f00*( z1^at + z2^at )^( 1/at-1 )*z2^( at-1 )*exp(z2)
  )
  
  f11 <- ( 
    f00*( ( z1^at + z2^at )^( 2/at-2 )*z1^( at-1 )*z2^( at-1 ) + ( z1^at + z2^at )^( 1/at-2 )*z1^( at-1 )*z2^( at-1 )*( at-1 ) )*exp(z1)*exp(z2)
  )
  
  llh <- sum(
    sapply(1:length(marginal_label), function(k){ 
      # if(is.na(get( paste0('f',marginal_label[k] ) )[k])){
      #   print( get( paste0('f',marginal_label[k] ) )[k] )
      #   print( yt )
      #   print( parameters )
      #   print( marginal_label[k] )
      #   print( k )
      # }
      return(
        # get( paste0('f',marginal_label[k] ) )[k] 
        log( get( paste0('f',marginal_label[k] ) )[k] )
      )} )
      )
  # if(is.na(llh)){
  #   browser()
  #   }
  # plot(at)
  return( list( llh = -llh, at = at, yt= yt ) )
}

fit_tvc <- function( marginal_variable, marginal_label){
  # 是否边际不在[0,1]之间
  for (i in length(marginal_variable)) {
    if( max(marginal_variable[[i]])>1 || min(marginal_variable[[i]])<0){
      print("边际不在[0,1]内")      
    }
  }
  
  Z <- list( -log(marginal_variable[[1]]), -log(marginal_variable[[2]]) )
  
  lf <-function(parameters){
    llh_tvc(Z, parameters, marginal_label)[[1]]
  }   
  par.start <- c(0.01, 0.5, 0.1)
  # ui.opt <- rbind(c(0, -1, -1),
  #                 c(0, 1,  0),
  #                 c(0, 0,  1))
  ui.opt <- rbind(c(0, -1, 0),
                  c(0, 1,  0),
                  c(0, 0,  1))
  ci.opt <- c(-0.999999, 0, 0)
  p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                            grad = NULL,
                            ui = ui.opt, ci = ci.opt,
                            hessian = FALSE)
  # p.e.nlminb <- optim( par.start, f = function(theta) { sum(lf(theta)) }, method = "BFGS")
  
  print(p.e.nlminb)
  p.e.nlminb$value <- -p.e.nlminb$value
  par <- p.e.nlminb$par
  
  # 计算时变参数
  at = llh_tvc(Z, par, marginal_label)[[2]]
  yt = llh_tvc(Z, par, marginal_label)[[3]]
  # 预测下一天的yt和at
  ut = calculate_yt( Z, par, marginal_label)[[2]]
  predict_yt = as.numeric( par%*%c( 1, tail( yt, 1), tail( ut, 1) ) )
  predict_at = 1 + exp(predict_yt)
  # 计算标准误和p值
  inv_hessian <- try({
    solve(-suppressWarnings(hessian(x = par, func = function (theta) {
      -sum(lf(theta))
    }, method.args=list(eps=1e-4, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))))
  }, silent = TRUE)


  if (class(inv_hessian)[1] == "try-error") {
    warning("Inverting the Hessian matrix failed. No robust standard errors calculated. Possible workaround: Multiply returns by 100.")
    rob.std.err <- NA
  } else {
    rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(func = lf, x = par)) %*% inv_hessian))
  }
  
  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         # std.err = rob.std.err,
         broom.mgarch = data.frame(
                                   # term = names(par),
                                   estimate = par
                                   , rob.std.err = rob.std.err
                                   , p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err))))
                                   ),
         at = at,
         yt = yt,
         ut = ut,
         predict_yt = predict_yt,
         predict_at = predict_at,
         cdf = exp( -(Z[[1]]^at+Z[[2]]^at)^(1/at) ),
         llh = p.e.nlminb$value,
         optim = p.e.nlminb
    )
  
  # Add class mfGARCH for employing generic functions
  class(output) <- "tvc"
  output
}

# gumbel的一阶偏导，对x1求偏导
pd_gumbel <- function( x1, x2, a){
  exp( -((-log(x1))^a+(-log(x2))^a)^(1/a) )*( (-log(x1))^a + (-log(x2))^a )^( 1/a-1 )*(-log(x1))^( a-1 )*x1^(-1)
}

gumbel <- function( x1, x2, a){
  exp( -((-log(x1))^a+(-log(x2))^a)^(1/a) )
}

# 以x1的p作为条件，parameters为x2即为目标变量的分布参数
CoVaR <- function( p, q, parameters, ret, threshold, a, family = 'gumbel', distribution = '2PoT', method){
  # VaR_x1 <- two_side_reverse_Pot_distribution( p, parameters, ret, threshold )
  if ( family == 'gumbel' ){
    # f <- function(theta) { (100*(pd_gumbel( p, theta, a)-q))^2 }
    f <- function(theta) { (100*( dduCopula( c( p, theta), gumbelCopula( param = a ) ) - q ) )^2 }
    
  }else if( family == 't' ){
    f <- function(theta) { (100*( dduCopula( c( p, theta), tCopula( param = a ) ) - q ) )^2 }
  }
  
  theta <-  optim(0.5, function(x) f(x), lower = 0, upper = 1)$par
  
  if( method == 1 ){
    if( distribution == '2PoT' ){
      PoT_CoVaR <- two_side_reverse_Pot_distribution( theta, parameters, ret, threshold )
    }else if( distribution == 'PoT' ){
      PoT_CoVaR <- one_side_reverse_Pot_distribution( theta, parameters, ret, threshold = threshold[[1]] )
    }
    result <- PoT_CoVaR
    if(PoT_CoVaR > 0.00001){
      result <- sort(ret)[ floor( (threshold[[1]])*length(ret) ) ] + PoT_CoVaR
    }else if( PoT_CoVaR < (-0.00001) ){
      result <- sort(ret)[ floor( (1-threshold[[1]])*length(ret) ) ] + PoT_CoVaR
    }  
  }else if ( method == 2 ){
    if( distribution == '2PoT' ){
      PoT_CoVaR <- two_side_reverse_residuals_distribution( theta, parameters, ret, threshold )
    }else if( distribution == 'PoT' ){
      PoT_CoVaR <- one_side_reverse_residuals_distribution( theta, parameters, ret, threshold = threshold[[1]] )
    }
    result <- PoT_CoVaR
  }
  
  return( result )
}

g_CoVaR <- function( p, q, parameters, ret, threshold, a, family = 'gumbel', distribution = '2PoT', method){
  if( p>0.5 && q>0.5 ){
    # f <- function(theta) { (100*( ( ( theta - gumbel( p, theta, a) )/p) - q ) )^2 }#上行风险参见On dependence consistency of CoVaR and some other systemic risk measures第7页
    if ( family == 'gumbel' ){
      f <- function(theta) { (100*( ( 1 - p - theta + pCopula( c( p, theta), gumbelCopula( param = a ) ) ) - (1-p)*(1-q) ) )^2 }
      
    }else if( family == 't' ){
      f <- function(theta) { (100*( ( 1 - p - theta + pCopula( c( p, theta), tCopula( param = a ) ) ) - (1-p)*(1-q) ) )^2 }
    }
  }else if( p<0.5 && q<0.5 ){
    # f <- function(theta) { (100*( (gumbel( p, theta, a)/p) - q ) )^2 }#下行风险

    if ( family == 'gumbel' ){
      f <- function(theta) { (100*( pCopula( c( p, theta), gumbelCopula( param = a ) ) -p*q ) )^2 }

    }else if( family == 't' ){
      f <- function(theta) { (100*( pCopula( c( p, theta), tCopula( param = a ) ) - p*q ) )^2 }
    }

  }
  # # VaR_x1 <- two_side_reverse_Pot_distribution( p, parameters, ret, threshold )
  

  
  theta <-  optim(0.5, function(x) f(x), lower = 0, upper = 1)$par
  if( method == 1 ){
    if( distribution == '2PoT' ){
      PoT_CoVaR <- two_side_reverse_Pot_distribution( theta, parameters, ret, threshold )
    }else if( distribution == 'PoT' ){
      PoT_CoVaR <- one_side_reverse_Pot_distribution( theta, parameters, ret, threshold = threshold[[1]] )
    }
    result <- PoT_CoVaR
    if(PoT_CoVaR > 0.00001){
      result <- sort(ret)[ floor( (threshold[[1]])*length(ret) ) ] + PoT_CoVaR
    }else if( PoT_CoVaR < (-0.00001) ){
      result <- sort(ret)[ floor( (1-threshold[[1]])*length(ret) ) ] + PoT_CoVaR
    }  
  }else if ( method == 2 ){
    if( distribution == '2PoT' ){
      PoT_CoVaR <- two_side_reverse_residuals_distribution( theta, parameters, ret, threshold )
    }else if( distribution == 'PoT' ){
      PoT_CoVaR <- one_side_reverse_residuals_distribution( theta, parameters, ret, threshold = threshold[[1]] )
    }
    result <- PoT_CoVaR
  }
  return( result )
}

# CoVaR回测检验
BacktestCoVaR<- function (data, VaR, alpha, x2_VaR_beta ) 
{
  # vY = data[1,]
  # 找到条件下的数据集
  if( alpha >0.5){
    coditional_set = which(data[2,] > x2_VaR_beta)
  }else {
    coditional_set = which(data[2,] < x2_VaR_beta)
  }
  # print(coditional_set)
  
  vY = data[1,][coditional_set]
  vVaR = VaR[coditional_set]
  dTau = alpha
  vY = as.numeric(vY)
  vVaR = as.numeric(vVaR)

  # index= which(vY[coditional_set]<VaR[coditional_set])
  if( alpha >0.5){
    # Hit = as.numeric( rep(1,length(vY)) )
    # Hit[coditional_set[which(vY[coditional_set]>vVaR[coditional_set])]] = 0
    Hit = as.numeric( rep(0,length(vY)) )
    Hit[vY<vVaR] = 1
  }
  else{
    # Hit = as.numeric( rep(0,length(vY)) )
    # Hit[coditional_set[which(vY[coditional_set]<vVaR[coditional_set])]] = 1
    Hit = as.numeric( rep(0,length(vY)) )
    Hit[vY<vVaR] = 1
  }
  # print(which(vY[coditional_set]>vVaR[coditional_set]))
  LRuc = Kupiec(Hit, dTau)
  # LRcc = Christoffersen(Hit, dTau)
  # DQ = DQOOStest(vY, vVaR, dTau, Lags)
  # AE = ActualOverExpected(Hit, dTau)
  # AD = AbsoluteDeviation(Hit, vY, vVaR)
  # Loss = QLoss(vY, vVaR, dTau)
  # lOut = list(LRuc = LRuc, LRcc = LRcc, AE = AE, AD = AD, 
  #             DQ = DQ, Loss = Loss)
  Loss = QLoss(vY, vVaR, dTau, Hit)
  lOut = list(LRuc = LRuc, Loss = Loss)
  return(lOut)
}
# 
Kupiec <- function(Hit, tau) {
  N = length(Hit)
  x = sum(Hit)
  rate = x/N
  test = -2 * log(((1 - tau)^(N - x) * tau^x)/((1 - rate)^(N - x) * rate^x))
  if (is.nan(test))
    test = -2 * ((N - x) * log(1 - tau) + x * log(tau) - (N - x) * log(1 - rate) - x * log(rate))
  # threshold = qchisq(alphaTest, df = 1)
  pvalue = 1 - pchisq(test, df = 1)
  
  LRpof = c(test, pvalue)
  names(LRpof) = c("Test", "Pvalue")
  return(LRpof)
}
# 
QLoss <- function(vY, vVaR, dTau, vHit) {
  vLoss = (vY - vVaR) * (dTau - vHit)
  dLoss = mean(vLoss)
  
  return(list(Loss = dLoss, LossSeries = vLoss))
}


# 基于多种分布的garch边际计算联合copula后的var
rolling_prediction <- function( i ) {
  if ( distribution == '2PoT' ){
    # 对每一个变量进行拟合
    for (j in 1:length(variable_list)) {
      
      DataFrame=variable_list[[j]][[1]][(i-window_length+1):i,]#一定要括号起来,不然是错误的
      numbers <- str_extract_all(variable_list[[j]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      # numbers[1] <- 0.97
      assign( paste0( 'fit_model_', j), Pot_fit_garch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                                      low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                                      , pot_k = TRUE, two_side = TRUE, gamma = TRUE ))
      ret <- get( paste0( 'fit_model_', j) )$df.fitted$return
      # 计算边际的Pot_y
      assign(
        paste0( 'fit_marginal_Pot_y_', j), calculate_marginal_label( get( paste0( 'fit_model_', j) ),
                                                                     threshold = c( threshold_0 = as.numeric(numbers[1]), threshold_1 = as.numeric(numbers[1]) ) )[["Pot_y"]]
      )
      # 计算Pot_y对应的区域,得到<负阈值，>阈值以及中间部分的标签
      assign(
        paste0( 'fit_marginal_Pot_y_label_', j), calculate_marginal_label( get( paste0( 'fit_model_', j) ),
                                                                           threshold = c( threshold_0 = as.numeric(numbers[1]), threshold_1 = as.numeric(numbers[1]) ) )[["label"]]
      )
      # 计算边际Pot_y的累积分布变量
      sigma <- get( paste0( 'fit_model_', j) )$df.fitted$tau * get( paste0( 'fit_model_', j) )$df.fitted$g
      e_0 <- get( paste0( 'fit_model_', j) )$e$e_0
      e_1 <- get( paste0( 'fit_model_', j) )$e$e_1
      
      if (( "k_0" %in% names(get( paste0( 'fit_model_', j) )$par) ) ){
        k_0 <- get( paste0( 'fit_model_', j) )$par[["k_0"]]
        k_1 <- get( paste0( 'fit_model_', j) )$par[["k_1"]]
      }else{
        k_0 = 1
        k_1 = 1
      }
      
      parameters <- list( k_0 = k_0, k_1 = k_1, e_0 = e_0, e_1 = e_1, sigma = sigma)
      
      assign(
        paste0( 'fit_marginal_Pot_y_cdf_', j), Pot_y_cdf( ret,parameters, threshold = c( threshold_0 = as.numeric(numbers[1]),
                                                                                         threshold_1 = as.numeric(numbers[1]) ) )
      )
      
      print(j)
    }
  }else if(  distribution == 'PoT' ){
    
    # 对每一个变量进行拟合
    for (j in 1:length(variable_list)) {
      
      DataFrame=variable_list[[j]][[1]][(i-window_length+1):i,]#一定要括号起来,不然是错误的
      numbers <- str_extract_all(variable_list[[j]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      # numbers[1] <- 0.97
      assign( paste0( 'fit_model_', j), Pot_fit_garch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                                      low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                                      , pot_k = TRUE, two_side = FALSE, gamma = TRUE ))
      ret <- get( paste0( 'fit_model_', j) )$df.fitted$return
      # 计算边际的Pot_y
      assign(
        paste0( 'fit_marginal_Pot_y_', j), calculate_marginal_label( get( paste0( 'fit_model_', j) ),
                                                                     threshold = as.numeric(numbers[1]) )[["Pot_y"]]  )   
      
      assign(
        paste0( 'fit_marginal_Pot_y_label_', j), calculate_marginal_label( get( paste0( 'fit_model_', j) ),
                                                                           threshold = c( threshold = as.numeric(numbers[1]) ) )[["label"]]
      )
      # 计算边际Pot_y的累积分布变量
      sigma <- get( paste0( 'fit_model_', j) )$df.fitted$tau * get( paste0( 'fit_model_', j) )$df.fitted$g
      e <- get( paste0( 'fit_model_', j) )$e
      
      if (( "k" %in% names(get( paste0( 'fit_model_', j) )$par) ) ){
        k <- get( paste0( 'fit_model_', j) )$par[["k"]]
      }else{
        k = 1
      }
      
      parameters <- list( k = k, e = e, sigma = sigma)
      assign(
        paste0( 'fit_marginal_Pot_y_cdf_', j), Pot_y_cdf( ret,parameters, threshold = c( threshold = as.numeric(numbers[1]) ) )
      )
      
      print(j)
    }
  }else if(  distribution == 'norm' ){
    
    for (j in 1:length(variable_list)){
      DataFrame=variable_list[[j]][[1]][(i-window_length+1):i,]#一定要括号起来,不然是错误的
      numbers <- str_extract_all(variable_list[[j]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      assign( paste0( 'best_fit_garch_', i), fit_mfgarch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])), low.freq = low.freq, K = K, distribution='norm') )
      assign( paste0( 'marginal_y_cdf_', i), na.omit(pnorm(get( paste0( 'best_fit_garch_', i) )$df.fitted$residuals ) ) )
      
    }
  }else if(  distribution == 't' ){
    
    for (j in 1:length(variable_list)){
      DataFrame=variable_list[[j]][[1]][(i-window_length+1):i,]#一定要括号起来,不然是错误的
      numbers <- str_extract_all(variable_list[[j]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      assign( paste0( 'best_fit_garch_', i), fit_mfgarch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])), low.freq = low.freq, K = K, distribution='t') )
      assign( paste0( 'marginal_y_cdf_', i), na.omit(pt(get( paste0( 'best_fit_garch_', i) )$df.fitted$residuals ) ) )
      
    }
  }
  
  if( (distribution == '2PoT')||(distribution == 'PoT') ){
    fit_marginal_Pot_y_label <- paste0( fit_marginal_Pot_y_label_1, fit_marginal_Pot_y_label_2)
    fit_copula <- fit_tvc(list( fit_marginal_Pot_y_cdf_1, fit_marginal_Pot_y_cdf_2), fit_marginal_Pot_y_label)
    
    # 随机生成copula
    # random_pair <- rCopula(10000, gumbelCopula( fit_copula$predict_at, dim = 2))
    random_pair <- rCopula(10000, gumbelCopula( fit_copula$predict_at, dim = 2))
  }else if(  (distribution == 'norm')||(distribution == 't') ){
    
    fit_copula <- BiCopSelect(marginal_y_cdf_1,marginal_y_cdf_2)
    # 随机生成copula
    random_pair <-  BiCopSim(10000, fit_copula)
  }
  
  
  if ( distribution == '2PoT' ){
    # 逆变换得到收益率
    for (w in 1:length(variable_list) ) {
      # 获取阈值、外生变量信息
      numbers <- str_extract_all(variable_list[[w]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      # assign( 
      #   paste0( 'numbers_', w), str_extract_all(variable_list[[w]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      # )
      # numbers <- get( paste0( 'numbers_', w) )
      # numbers[1] <- 0.97
      # 初始化随机收益率数组
      assign( paste0('random_xt_',w), rep( NA, nrow(random_pair) ) )
      # 初始化(预测下一天)边际上的参数
      fit_model <- get( paste0( 'fit_model_', w) )
      e_0 <- calculate_shape_value( fit_model$par[['Pot_p1_0']], fit_model$par[['Pot_p2_0']], fit_model$par[['Pot_p3_0']]
                                    , c(fit_model$df.fitted$return,0) )# 多一个长度是为了预测t+1时刻的尾部参数
      e_1 <- calculate_shape_value( fit_model$par[['Pot_p1_1']], fit_model$par[['Pot_p2_1']], fit_model$par[['Pot_p3_1']]
                                    , c(fit_model$df.fitted$return,0) )
      if (( "k_0" %in% names(fit_model$par) ) ){
        k_0 <- fit_model$par[["k_0"]]
        k_1 <- fit_model$par[["k_1"]]
      }else{
        k_0 = 1
        k_1 = 1
      }
      # 预测下一天的参数以构建下一天的收益率
      # parameters <- c( sigma = predict( fit_model, c(1) ), e_0 = e_0[length(e_0)], e_1 = e_1[length(e_1)], k_0 = k_0, k_1 = k_1 )
      assign( paste0( 'parameters_', w), c( sigma = predict( fit_model, c(1) ), e_0 = e_0[length(e_0)], e_1 = e_1[length(e_1)], k_0 = k_0, k_1 = k_1 ) )
   
      if( method == 1 ){
        assign( paste0('random_Pot_',w),sapply( 1:nrow(random_pair), function(q){
          two_side_reverse_Pot_distribution( random_pair[q,w]
                                             , get( paste0( 'parameters_', w) ), fit_model$df.fitted$return
                                             , threshold = c( threshold_0 = as.numeric( numbers[1] )
                                                              , threshold_1 = as.numeric( numbers[1] ) ) )
        } ) )
        
        assign( paste0('random_xt_',w),
                ifelse(get(paste0('random_Pot_',w)) > 0.001,
                       sort(fit_model$df.fitted$return)[ floor( ( as.numeric( numbers[1] ) )*length(fit_model$df.fitted$return) ) ] + get(paste0('random_Pot_',w)),
                       ifelse( get( paste0('random_Pot_',w)) < (-0.001),
                               sort(fit_model$df.fitted$return)[ floor( ( 1-as.numeric( numbers[1] ) )*length(fit_model$df.fitted$return) ) ] + get(paste0('random_Pot_',w)),
                               get(paste0('random_Pot_',w))
                       )
                )
        )
      }else if( method == 2 ){
        # # 逆分布为收益率
        assign( paste0('random_xt_',w),sapply( 1:nrow(random_pair), function(j){
                                                two_side_reverse_residuals_distribution(random_pair[j,w]
                                                , get( paste0( 'parameters_', w) ), fit_model$df.fitted$return
                                                , threshold = as.numeric( numbers[1] ))}) )
      }
 
    }
  }else if( distribution == 'PoT' ){
    # 逆变换得到收益率
    for (w in 1:length(variable_list) ) {
      # 获取阈值、外生变量信息
      numbers <- str_extract_all(variable_list[[w]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      # 初始化随机收益率数组
      assign( paste0('random_xt_',w), rep( NA, nrow(random_pair) ) )
      # 初始化(预测下一天)边际上的参数
      fit_model <- get( paste0( 'fit_model_', w) )
      e <- calculate_shape_value( fit_model$par[['Pot_p1']], fit_model$par[['Pot_p2']], fit_model$par[['Pot_p3']]
                                  , c(fit_model$df.fitted$return,0) )# 多一个长度是为了预测t+1时刻的尾部参数
      if (( "k" %in% names(fit_model$par) ) ){
        k <- fit_model$par[["k"]]
      }else{
        k = 1
      }
      # 预测下一天的参数以构建下一天的收益率
      assign( paste0( 'parameters_', w), c( sigma = predict( fit_model, c(1) ), e = e[length(e)], k = k ) )
      # # 逆分布为收益率
      if( method == 1 ){
        assign( paste0('random_Pot_',w),sapply( 1:nrow(random_pair), function(q){
          one_side_reverse_Pot_distribution( random_pair[q,w]
                                             , get( paste0( 'parameters_', w) ), fit_model$df.fitted$return
                                             , threshold = c( threshold = as.numeric(numbers[1]) ) )
        } ) )
        
        assign( paste0('random_xt_',w),
                ifelse(get(paste0('random_Pot_',w)) > 0.001,
                       sort(fit_model$df.fitted$return)[ floor( ( as.numeric( numbers[1] ) )*length(fit_model$df.fitted$return) ) ] + get(paste0('random_Pot_',w)),
                       get(paste0('random_Pot_',w))
                )
        )
      }else if( method == 2 ){
        assign( paste0('random_xt_',w),sapply( 1:nrow(random_pair), function(j){
          one_side_reverse_residuals_distribution(random_pair[j,w]
                                            , get( paste0( 'parameters_', w) ), fit_model$df.fitted$return
                                            , threshold = c( threshold = as.numeric(numbers[1]) ) )}) )
      }
    }
  }else if( distribution == 'norm' ){
    
    # 逆变换得到收益率
    for (w in 1:length(variable_list) ) {
      
      fit_model <- get( paste0( 'fit_model_', w) )
      assign( paste0('random_et_',w ), qnorm(random_pair[,w]) )
      assign( paste0('random_xt_',w ), get( paste0('random_et_',w) )*sqrt(predict( fit_model, c(1) ) ) )
      
      
    }
    
  }else if( distribution == 't' ){
    
    # 逆变换得到收益率
    for (w in 1:length(variable_list) ) {
      
      fit_model <- get( paste0( 'fit_model_', w) )
      assign( paste0('random_et_',w ), qt(random_pair[,w]) )
      assign( paste0('random_xt_',w ), get( paste0('random_et_',w) )*sqrt(predict( fit_model, c(1) , method = 'unnormal') ) )
      
      
    }
    
  }
  
  random_xt <- 0.5*random_xt_1 + 0.5*random_xt_2
  
  for(s in 1:length(VaR_level) ){
    VaR_set$total_VaR_set_up[s] <- sort(random_xt)[ floor( VaR_level[s]*length(random_xt) ) ]# 最终计算的var是总的收益率的样本分位数
    VaR_set$total_VaR_set_down[s]<- sort(random_xt)[ floor( ( 1 - VaR_level[s] )*length(random_xt) ) ]# 最终计算的var是总的收益率的样本分位数
    
    VaR_set$X1_VaR_set_up[s] <- sort(random_xt_1)[ floor( VaR_level[s]*length(random_xt_1) ) ]
    VaR_set$X1_VaR_set_down[s] <- sort(random_xt_1)[ floor( ( 1 - VaR_level[s] )*length(random_xt_1) ) ]
    VaR_set$X1_CoVaR_set_up[s] <- CoVaR( VaR_level[s], VaR_level[s],parameters_1, fit_model_1$df.fitted$return ,
                                           threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                          threshold_1 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                             a = fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X1_CoVaR_set_down[s] <- CoVaR( (1-VaR_level[s]), (1-VaR_level[s]),parameters_1, fit_model_1$df.fitted$return ,
                                             threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                            threshold_1 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                             a = fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X1_g_CoVaR_set_up[s] <- g_CoVaR( VaR_level[s], VaR_level[s],parameters_1, fit_model_1$df.fitted$return ,
                                               threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                              threshold_1 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                             a = fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X1_g_CoVaR_set_down[s] <- g_CoVaR( (1-VaR_level[s]), (1-VaR_level[s]),parameters_1, fit_model_1$df.fitted$return ,
                                                 threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                                threshold_1 = as.numeric(str_extract_all(variable_list[[1]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                                a = fit_copula$predict_at, distribution = distribution, method = method)
    
    
    VaR_set$X2_VaR_set_up[s] <- sort(random_xt_2)[ floor( VaR_level[s]*length(random_xt_2) ) ]
    VaR_set$X2_VaR_set_down[s] <- sort(random_xt_2)[ floor( ( 1 - VaR_level[s] )*length(random_xt_2) ) ]
    VaR_set$X2_CoVaR_set_up[s] <- CoVaR( VaR_level[s], VaR_level[s],parameters_2, fit_model_2$df.fitted$return ,
                                           threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                          threshold_1 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                         a = fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X2_CoVaR_set_down[s] <- CoVaR( (1-VaR_level[s]), (1-VaR_level[s]),parameters_2, fit_model_2$df.fitted$return ,
                                             threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                            threshold_1 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                         a = fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X2_g_CoVaR_set_up[s] <- g_CoVaR( VaR_level[s], VaR_level[s],parameters_2, fit_model_2$df.fitted$return ,
                                               threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                              threshold_1 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                         a =  fit_copula$predict_at, distribution = distribution, method = method)
    VaR_set$X2_g_CoVaR_set_down[s] <- g_CoVaR( (1-VaR_level[s]), (1-VaR_level[s]),parameters_2, fit_model_2$df.fitted$return ,
                                                 threshold = c( threshold_0 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]),
                                                                threshold_1 = as.numeric(str_extract_all(variable_list[[2]][[2]], "\\d+\\.?\\d*")[[1]][[1]]) ),
                                         a = fit_copula$predict_at, distribution = distribution, method = method)
    # Loop through VaR_set and check conditions to assign values
    for (ss in 1:length(VaR_set)) {
      if (grepl("total", names(VaR_set)[ss]) & grepl("up", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1) {
        VaR_set[[ss]][[s]] <- sort(fit_model_1$df.fitted$return + fit_model_2$df.fitted$return)[ floor( as.numeric( numbers[1] )*length(fit_model_1$df.fitted$return)) ]
      }else if (grepl("total", names(VaR_set)[ss]) & grepl("down", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1 ) {
        VaR_set[[ss]][[s]] <- sort(fit_model_1$df.fitted$return + fit_model_2$df.fitted$return)[ floor( ( 1-as.numeric( numbers[1] ) )*length(fit_model_1$df.fitted$return)) ]
      }else if (grepl("X1", names(VaR_set)[ss]) & grepl("up", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1 ) {
        VaR_set[[ss]][[s]] <- sort(fit_model_1$df.fitted$return)[ floor( as.numeric( numbers[1] )*length(fit_model_1$df.fitted$return) ) ]
      }else if (grepl("X1", names(VaR_set)[ss]) & grepl("down", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1 ) {
        VaR_set[[ss]][[s]] <- sort(fit_model_1$df.fitted$return)[ floor( ( 1-as.numeric( numbers[1] ) )*length(fit_model_1$df.fitted$return) ) ]
      }else if (grepl("X2", names(VaR_set)[ss]) & grepl("up", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1 ) {
        VaR_set[[ss]][[s]] <- sort(fit_model_2$df.fitted$return)[ floor( as.numeric( numbers[1] )*length(fit_model_2$df.fitted$return) ) ]
      }else if (grepl("X2", names(VaR_set)[ss]) & grepl("down", names(VaR_set)[ss]) & abs(VaR_set[[ss]][[s]]) < 0.1 ) {
        VaR_set[[ss]][[s]] <- sort(fit_model_2$df.fitted$return)[ floor( ( 1-as.numeric( numbers[1] ) )*length(fit_model_2$df.fitted$return) ) ]
      }
    }
     
  }
  

  
  return( list( VaR_set = VaR_set, parameters = c( sigma_1 = predict( fit_model_1, c(1) ), sigma_2 = predict( fit_model_2, c(1) ), at = fit_copula$predict_at ) ) )
}
# 回测检验

BacktestVaR_result <- function(total_return, total_VaR_set_up, total_VaR_set_down) {
  # 适用BacktestVaR进行回测
  results_up <- sapply(1:5, function(i) BacktestVaR(total_return, total_VaR_set_up[,i][!is.na(total_VaR_set_up[,i])], c(0.9, 0.95, 0.975, 0.99, 0.995)[i])$LRuc[[2]])
  results_down <- sapply(1:5, function(i) BacktestVaR(-(total_return), -(total_VaR_set_down[,i][!is.na(total_VaR_set_down[,i])]), c(0.9, 0.95, 0.975, 0.99, 0.995)[i])$LRuc[[2]])
                         

  return(
    data.frame(
      level = as.numeric(c(0.9, 0.95, 0.975, 0.99, 0.995)),results_up = results_up,
      # ratio_up = round(sapply(1:5, function(i) sum((total_return - total_VaR_set_up[,i]) < 0) / length(total_return)), 6),
      results_down = results_down
      # , ratio_down = round(sapply(1:5, function(i) sum((total_return - total_VaR_set_down[,i]) < 0) / length(total_return)), 6)            
              )
       )
}

BacktestCoVaR_result <- function(return_1, return_2, X1_CoVaR_set_up, X2_VaR_set_up, X1_CoVaR_set_down, X2_VaR_set_down) {
  results_up <- sapply(1:5, function(i) BacktestCoVaR(rbind(return_1, return_2), X1_CoVaR_set_up[,i][!is.na(X1_CoVaR_set_up[,i])], c(0.9, 0.95, 0.975, 0.99, 0.995)[i], X2_VaR_set_up[,i][!is.na(X2_VaR_set_up[,i])])$LRuc[[2]])
  results_down <- sapply(1:5, function(i) BacktestCoVaR(rbind(return_1, return_2), X1_CoVaR_set_down[,i][!is.na(X1_CoVaR_set_down[,i])], c(0.1, 0.05, 0.025, 0.01, 0.005)[i], X2_VaR_set_down[,i][!is.na(X2_VaR_set_down[,i])])$LRuc[[2]])

  return(data.frame(level = as.numeric(c(0.9, 0.95, 0.975, 0.99, 0.995)),results_up = results_up,results_down = results_down))
}

# 绘图
create_plot <- function(X1_VaR_set_up, X1_VaR_set_down, total_return, method = '1',distribution = '2PoT',data_tail = 'up', VaR) {
  data <- data.frame( level_995 = X1_VaR_set_up[,5][!is.na(X1_VaR_set_up[,5])],
                      level_99 = X1_VaR_set_up[,4][!is.na(X1_VaR_set_up[,4])],
                      level_975 = X1_VaR_set_up[,3][!is.na(X1_VaR_set_up[,3])],
                      level_95 = X1_VaR_set_up[,2][!is.na(X1_VaR_set_up[,2])],
                      level_90 = X1_VaR_set_up[,1][!is.na(X1_VaR_set_up[,1])],
                      level_005 = X1_VaR_set_down[,5][!is.na(X1_VaR_set_down[,5])],
                      level_01 = X1_VaR_set_down[,4][!is.na(X1_VaR_set_down[,4])],
                      level_025 = X1_VaR_set_down[,3][!is.na(X1_VaR_set_down[,3])],
                      level_05 = X1_VaR_set_down[,2][!is.na(X1_VaR_set_down[,2])],
                      level_1 = X1_VaR_set_down[,1][!is.na(X1_VaR_set_down[,1])],
                      real = total_return )
  if( distribution == '2PoT' ){
    p <- ggplot(data, aes(x = 1:nrow(data))) +
      geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
      # geom_line(aes(y = VaR_up_995, color = 'VaR995'), size = 0.7) +
      # geom_line(aes(y = VaR_up_99, color = 'VaR99'), size = 0.7) +
      geom_line(aes(y = level_975, color = 'level_975'), size = 0.7) +
      geom_line(aes(y = level_95, color = 'level_95'), size = 0.7) +
      # geom_line(aes(y = VaR_up_90, color = 'VaR90'), size = 0.7) +
      # geom_line(aes(y = VaR_down_995, color = 'VaR995'), size = 0.7) +
      # geom_line(aes(y = VaR_down_99, color = 'VaR99'), size = 0.7) +
      geom_line(aes(y = level_025, color = 'level_025'), size = 0.7) +
      geom_line(aes(y = level_05, color = 'level_05'), size = 0.7) +
      # geom_line(aes(y = VaR_down_90, color = 'VaR90'), size = 0.7) +
      scale_color_manual(values = c('real' = 'red', 'level_975' = 'blue', 'level_95' = 'purple', 'level_025' = 'blue', 'level_05' = 'purple') )+
      # labs(title = paste0( ifelse(method=='2','传统','改进'),'方法','下',distribution,'的',VaR), x = "时间", y = "收益率")
      labs(title = paste0( distribution,' ',VaR), x = "time", y = "return")
    
    p + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) # 背景白底加标题居中  
  }else if( ( distribution == 'PoT' )&&( data_tail == 'up' ) ){
    p <- ggplot(data, aes(x = 1:nrow(data))) +
      geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
      geom_line(aes(y = level_975, color = 'level_975'), size = 0.7) +
      geom_line(aes(y = level_95, color = 'level_95'), size = 0.7) +
      scale_color_manual(values = c('real' = 'red', 'level_975' = 'blue', 'level_95' = 'purple'))+
      # labs(title = paste0( ifelse(method=='2','传统','改进'),'方法','下',distribution,'_',data_tail,'的',VaR), x = "时间", y = "收益率")
      labs(title = paste0( distribution,'_',data_tail,' ',VaR), x = "time", y = "return")
    
    p + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) # 背景白底加标题居中 
  }else if( ( distribution == 'PoT' )&&( data_tail == 'down' ) ){
    data <- -(data)
    p <- ggplot(data, aes(x = 1:nrow(data))) +
      geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
      geom_line(aes(y = level_975, color = 'level_025'), size = 0.7) +
      geom_line(aes(y = level_95, color = 'level_05'), size = 0.7) +
      scale_color_manual(values = c('real' = 'red', 'level_025' = 'blue', 'level_05' = 'purple'))+
      # labs(title = paste0( ifelse(method=='2','传统','改进'),'方法','下',distribution,'_',data_tail,'的',VaR), x = "时间", y = "收益率")
      labs(title = paste0( distribution,'_',data_tail,' ',VaR), x = "time", y = "return")
    
    p + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) # 背景白底加标题居中 
    
  }
  
}
