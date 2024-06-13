#如果初始的e值太小会导致尾部概率大于1，因此  e[0] = exp(b1)是不对的;
cppFunction('NumericVector calculate_shape_value(double b1, double b2, double b3, NumericVector y) {
  int n = y.size();
  NumericVector e(n, NA_REAL);
  e[0] = exp(b1);
  
  for (int i = 1; i < n; i++) {
    e[i] = exp(b1 + b2 * log(e[i-1]) + b3 * exp(-fabs(y[i-1])));
  }
  
  e[0] = mean(e);
  
  for (int i = 1; i < n; i++) {
    e[i] = exp(b1 + b2 * log(e[i-1]) + b3 * exp(-fabs(y[i-1])));
  }
  
  return e;
}')

cppFunction(
  'double sum_tau_fcts(int i, double m, double theta, NumericVector phivar, NumericVector covariate, int K) {

  double exponential = m;

  for (int j = 1; j <= K; j++) {
    exponential += theta * phivar[j-1] * covariate[i - 1 - j];
  }

  return exponential;
}')

cppFunction(
  'NumericVector sum_tau(double m, double theta, NumericVector phivar, NumericVector covariate, int K) {

  int n = covariate.size() - K;
  NumericVector exponential(n);

  for (int i = 1; i <= n; i++) {
    exponential[i-1] = m;
    for (int j = 1; j <= K; j++) {
      exponential[i-1] += theta * phivar[j-1] * covariate[K + i - 1 - j];
    }
  }
  return exponential;
}')

cppFunction('NumericVector calculate_g(double omega, double alpha, double beta, double gamma, NumericVector returns, double g0) {
  int n = returns.size();

  NumericVector g(n);
  g[0] = g0;


  for (int i = 1; i < n; i++) {
    if (returns[i-1] >= 0) {
      g[i] = omega + alpha * pow(returns[i-1], 2) + beta * g[i-1];
    } else {
      g[i] = omega + alpha * pow(returns[i-1], 2) + gamma * pow(returns[i-1], 2) + beta * g[i-1];
    }
  }
  return g;
}')


#' @keywords internal
forecast_garch <- function(omega, alpha, beta, gamma, g, ret, steps.ahead) {
  omega / (1 - alpha - gamma/2 - beta) + (alpha + beta + gamma/2)^(steps.ahead - 1) * (omega + (alpha + gamma/2 * as.numeric(ret < 0)) * ret^2 + beta * g - omega / (1 - alpha - gamma/2 - beta))
}


# Calculates the phi weighting scheme
#' @keywords internal
calculate_phi <- function(w1, w2, K) {
  weights <- sapply(c(1:K),
                    FUN = function(j) (j / (K + 1))^(w1 - 1) * (1 - j / (K + 1))^(w2 - 1))
  weights <- weights/sum(weights)
  weights
}


#' @keywords internal
calculate_tau <- function(covariate, w1, w2, theta, m, K) { # used for simulation
  phi_var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), covariate)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi_var, covariate = covariate, K = K)))
  tau
}

#' @keywords internal
calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K,
                             x.two = NULL, K.two = NULL, theta.two = NULL,
                             low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))
  
  result <- merge(df, cbind(unique(df[low.freq]), tau), by = low.freq)
  
  if (!is.null(x.two)) {
    phi.var.two <- calculate_phi(w1.two, w2.two, K.two)
    covariate.two <- c(rep(NA, times = K.two), x.two)
    tau.two <- c(rep(NA, times = K.two),
                 exp(sum_tau(m = 0, theta = theta.two, phivar = phi.var.two,
                             covariate = x.two, K = K.two)))
    result <- merge(result, cbind(unique(df[low.freq.two]), tau.two), by = low.freq.two)
    
    result$tau.one <- result$tau # store tau component due to first covariate
    result$tau <- result$tau.one * result$tau.two # generate joint tau component
  }
  
  result
}


# # 为什么calculate_g函数没有文献中的公式tao，是因为输入参数return中已经除以根号tao
# calculate_g <- function(omega, alpha, beta, gamma, returns, g0) {
#   n <- length(returns)
#   g <- numeric(n)
#   g[1] <- g0
#   
#   for (i in 2:n) {
#     if (returns[i-1] >= 0) {
#       g[i] <- omega + alpha * returns[i-1]^2 + beta * g[i-1]
#     } else {
#       g[i] <- omega + alpha * returns[i-1]^2 + gamma * returns[i-1]^2 + beta * g[i-1]
#     }
#   }
#   
#   return(g)
# }


# df_test <- df_mfgarch %>% filter(is.na(vix) == FALSE)
# llh_mf(df_test, y = df_test$return, x = unlist(unique(df_test[c("date", "vix")])["vix"]), K = 3, mu = 0, m = 0, theta = 0.1, low.freq = "date", omega = 0.01, alpha = 0.06, beta = 0.9, gamma = 0, g_zero = 1, x.two = unlist(unique(df_test[c("year_month", "dindpro")])["dindpro"]), K.two = 12, theta.two = 0.1, low.freq.two = "year_month", w1.two = 1, w2.two = 1)
#' @keywords internal
llh_mf <-
  function(df, x, y, low.freq, mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (!is.null(x.two)) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    
    ret <- y
    ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood # is.na返回逻辑值trueorflase，1和0，which.min返回最小值
    tau <- tau[which.min(is.na(tau)):length(tau)]  # 因此0的最小下表就是tau有值的第一个位置的索引
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                     returns = ((ret - mu)/sqrt(tau)), g0 = g_zero)
    
    if (sum(g <= 0) > 0) {
      #rep(NA, times = length(y))
      #stop("g_t seems to be negative for at least one point in time?")
      rep(NA, times = length(g))
    } else {
      1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
    }
  }



#那么似然函数将变为
llh_mf_Pot <-
  function(df, x, y,Pot_y, low.freq,mu=NULL,
           k = 1,
           omega, alpha, beta, gamma,
           m, theta, 
           # positive_threshold_value,negative_threshold_value,
           threshold_value,
           Pot_p1, Pot_p2, Pot_p3, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (!is.null(x.two)) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    ###
    # k=1
    ret <- y#-mu
    ret <- ret[which.min(is.na(tau)):length(tau)]  # lags can't be used for likelihood # is.na返回逻辑值trueorflase，1和0，which.min返回最小值
    Pot_y <- Pot_y[which.min(is.na(tau)):length(tau)]#其实严格意义上Pot_y应该是在llh_mf_pot中计算的，
                                                     #因为研究的样本长度应该是舍弃k*low.freq后的长度，需要计算两边的阈值再计算Pot_y
    tau <- tau[which.min(is.na(tau)):length(tau)]  # 因此0的最小下标就是tau有值的第一个位置的索引
    yt = sign(Pot_y)
    
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma, returns = ((ret)/sqrt(tau)),# g0 = g_zero
                     g0 = g_zero)
    sigma1<-sqrt(tau*g)
    if(length(threshold_value)>1){
      e_0<-calculate_shape_value(b1 = Pot_p1[1],b2 = Pot_p2[1],b3 = Pot_p3[1],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      e_1<-calculate_shape_value(b1 = Pot_p1[2],b2 = Pot_p2[2],b3 = Pot_p3[2],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      # e<-calculate_shape_value(b1 = Pot_p1,b2 = Pot_p2,b3 = Pot_p3,(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      # penalty <- ifelse( ( 1-k[1]*(threshold_value[1]/sigma1)^(-e)-k[2]*(threshold_value[2]/sigma1)^(-e) )>0
      #                    , log( 1-k[1]*(threshold_value[1]/sigma1)^(-e)-k[2]*(threshold_value[2]/sigma1)^(-e) ),-1e10)
      penalty0 <- ifelse( ( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) )>0
                         , log( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) ),-1e10)
      penalty1 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                          , ( log(k[1])+log(e_0)+e_0*log(sigma1)-(e_0+1)*log(threshold_value[1]+Pot_y) )
                          , -1e6 )
      penalty2 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                          , ( log(k[2])+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]+Pot_y) )
                          , -1e6 )# 注意：当ret<0时理论上密度函数的表达式是
                                  # ( log(k[2])+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]-Pot_y) )
                                  # ,但是代码中Pot_y做了镜像处理，所以和ret>0是一样的表达式
      -(((yt-1)/-1)*penalty0
        +yt*(ret>0)*(penalty1)
        +yt*(ret<0)*(penalty2))#注意区分小k和大K
      # -(((yt-1)/-1)*penalty0
      # +yt*(ret>0)*(log(k[1])+log(e_0)+e_0*log(sigma1)-(e_0+1)*log(threshold_value[1]+Pot_y))
      # +yt*(ret<0)*(log(k[2])+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]+Pot_y)))#注意区分小k和大K

      # -(((yt-1)/-1)*penalty
      #   +yt*(ret>0)*(log(k[1])+log(e)+e*log(sigma1)-(e+1)*log(threshold_value[1]+Pot_y))
      #   +yt*(ret<0)*(log(k[2])+log(e)+e*log(sigma1)-(e+1)*log(threshold_value[2]+Pot_y)))#注意区分小k和大K

      # e_0<-calculate_shape_value(b1 = Pot_p1[1],b2 = Pot_p2[1],b3 = Pot_p3[1],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      # e_1<-calculate_shape_value(b1 = Pot_p1[2],b2 = Pot_p2[2],b3 = Pot_p3[2],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      # penalty0 <- ifelse( ( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) )>0
      #                     , log( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) ),-1e6 )
      # penalty1 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
      #                     , ( log(k[1])+log(e_0)+e_0*log(sigma1)-(e_0+1)*log(threshold_value[1]+Pot_y) )
      #                     , 0 )
      # penalty2 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
      #                     , ( log(k[2])+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]+Pot_y) )
      #                     , 0 )

    }else{
      e<-calculate_shape_value(b1 = Pot_p1,b2 = Pot_p2,b3 = Pot_p3,(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
      penalty0<- ifelse( ( 1-k*(threshold_value/sigma1)^(-e) )>0, log( 1-k*(threshold_value/sigma1)^(-e) ), -1e6 )
      penalty1 <- ifelse( 
        # ( ( k*( (threshold_value) /sigma1)^(-e) ) + k*( (threshold_value) /sigma1)^(-e)  )<1
        ( ( k*( (threshold_value) /sigma1)^(-e) ) )<1, ( log(k)+log(e)+e*log(sigma1)-(e+1)*log(threshold_value+Pot_y) )
                          , -1e6 )
      
      # print(penalty )
      if (sum(g <= 0) > 0) {
        #rep(NA, times = length(y))
        #stop("g_t seems to be negative for at least one point in time?")
        rep(NA, times = length(g))
      } else {
        -(((yt-1)/-1)*penalty0+yt*penalty1)#注意区分小k和大K
        # -(((yt-1)/-1)*log(abs(1-k*(threshold_value/sigma1)^(-e)))+yt*(log(k)+log(e)+e*log(sigma1)-(e+1)*log(threshold_value+Pot_y)))#注意区分小k和大K
      }
    }
  }


#老的极值理论
old_llh_mf_Pot <-
  function(df, x, y, low.freq,#, mu
           omega, alpha, beta, gamma,
           m, theta, threshold_value, Pot_p1, Pot_p2, Pot_p3, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (!is.null(x.two)) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    ###
    ret <- y#-mu
    #编写Pot_y取大于0的值
    threshold_value=sort(ret)[floor(length(ret)*threshold_value)]
    Pot_y = (ret-threshold_value)
    Pot_y = ifelse(Pot_y>0,Pot_y,0)
    
    ret <- ret[which.min(is.na(tau)):length(tau)]  # lags can't be used for likelihood # is.na返回逻辑值trueorflase，1和0，which.min返回最小值
    Pot_y = Pot_y[which.min(is.na(tau)):length(tau)]
    tau <- tau[which.min(is.na(tau)):length(tau)]  # 因此0的最小下表就是tau有值的第一个位置的索引
    yt = sign(Pot_y)
    
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = 0, returns = ((ret)/sqrt(tau)),# g0 = g_zero
                     g0 = g_zero)
    sigma1<-sqrt(tau*g)
    e<-calculate_shape_value(b1 = Pot_p1,b2 = Pot_p2,b3 = Pot_p3,y=(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
    
    if (sum(g <= 0) > 0) {
      #rep(NA, times = length(y))
      #stop("g_t seems to be negative for at least one point in time?")
      rep(NA, times = length(g))
    } else {
      -(((yt-1)/-1)*log(1-(threshold_value+1)^(-e))+yt*(log(sigma1)*e+log(e)-e*log(1+threshold_value)-(e+1)*log(sigma1+Pot_y)))#注意区分小k和大K
    }
  }



#' @keywords internal
llh_mf_t <-
  function(df, x, y, low.freq, V,mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (!is.null(x.two)) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    
    ret <- y-mu
    ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood # is.na返回逻辑值trueorflase，1和0，which.min返回最小值
    tau <- tau[which.min(is.na(tau)):length(tau)]  # 因此0的最小下表就是tau有值的第一个位置的索引
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                     returns = ((ret)/sqrt(tau)), g0 = g_zero)
    sigma1 <- g*tau
    if (sum(g <= 0) > 0) {
      #rep(NA, times = length(y))
      #stop("g_t seems to be negative for at least one point in time?")
      rep(NA, times = length(g))
    } else {
      -(log(gamma((V+1)/2))-log(gamma(V/2))-1/2*log((V-2)*pi)
      -((V+1)/2*log(1+(ret)^2/((V-2)*sigma1))+1/2*log(sigma1)))
      }
  }
#' @keywords internal
llh_simple <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu) / sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
}


#定义plot函数
#' @importFrom graphics lines
#' @export
plot.mfGARCH <- function(x, ...) {
  if (class(x) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }
  
  if (x$K == 0) {
    plot(x = x$df.fitted["date"], y = (x$df.fitted$g),
         type = "l",
         xlab = colnames(x$df.fitted)[3], ylab = "vol",
         main = "sqrt(g)", sub = "")
  } else {
    
    df_plot <- aggregate(x$df.fitted, by = list(x$df.fitted[, 3]), FUN = mean)
    #按周（也就是x$df.fitted[, 3]）分组，对每一周的所有数据去平局值
    if ( ('e' %in% colnames(df_plot) )&&('kk' %in% names(x$par) ) ) {
      y_max = max( max(((x$par[['k']])^((df_plot$e)^(-1)))*(df_plot$g*df_plot$tau),na.rm = TRUE) , max((df_plot$tau),na.rm = TRUE) )
      y_min = min( min(((x$par[['k']])^((df_plot$e)^(-1)))*(df_plot$g*df_plot$tau),na.rm = TRUE) , min((df_plot$tau),na.rm = TRUE) )
      par(mfrow = c(4, 1))
      #第一张图片
      plot(x = df_plot[, 1], y = ((x$par[['k']])^((df_plot$e)*(-1)))*(df_plot$g*df_plot$tau),
           # plot(x = df_plot[, 1], y = sqrt(df_plot$g*df_plot$tau),
           type = "l",
           xlab = colnames(x$df.fitted)[3], 
           # ylab = "vol",
           main = "sqrt(tau * g) and sqrt(tau) in red", sub = "",
           ylim = c(y_min, y_max)
      )
      # print(((x$par[['k']])^((df_plot$e)^(-1)))*sqrt(df_plot$g*df_plot$tau))
      # lines(x = df_plot[, 1],
      #       y = sqrt(df_plot$g*df_plot$tau),
      #       col = "blue")
      lines(x = df_plot[, 1],
            y = (df_plot$tau),
            col = "red")
      #第二张图片
      y_max = max( max((df_plot$g*df_plot$tau),na.rm = TRUE) , max((df_plot$tau),na.rm = TRUE) )
      y_min = min( min((df_plot$g*df_plot$tau),na.rm = TRUE) , min((df_plot$tau),na.rm = TRUE) )
      # plot(x = df_plot[, 1], y = ((x$par[['k']])^((df_plot$e)^(-1)))*sqrt(df_plot$g*df_plot$tau),
      plot(x = df_plot[, 1], y = (df_plot$g*df_plot$tau),
           type = "l",
           xlab = colnames(x$df.fitted)[3], 
           # ylab = "vol",
           # main = "sqrt(tau * g) and sqrt(tau) in red", sub = "",
           ylim = c(y_min, y_max)
      )
      # print(((x$par[['k']])^((df_plot$e)^(-1)))*sqrt(df_plot$g*df_plot$tau))
      lines(x = df_plot[, 1],
            y = (df_plot$tau),
            col = "red")
      #第三张图片
      plot(x = df_plot[, 1], y = ((df_plot$g))*((x$par[['k']])^((df_plot$e)^(-1))),
           type = "l")
      #第四张图片
      plot(x = df_plot[, 1], y = ((df_plot$g)),
           type = "l")
    } else {
      # par(mfrow = c(1, 1))
      # y_max = max( max(sqrt(df_plot$g*df_plot$tau),na.rm = TRUE) , max(sqrt(df_plot$tau),na.rm = TRUE) )
      # y_min = min( min(sqrt(df_plot$g*df_plot$tau),na.rm = TRUE) , min(sqrt(df_plot$tau),na.rm = TRUE) )
      # # plot(x = df_plot[, 1], y = ((x$par[['k']])^((df_plot$e)^(-1)))*sqrt(df_plot$g*df_plot$tau),
      # plot(x = df_plot[, 1], y = sqrt(df_plot$g*df_plot$tau),
      #      type = "l",
      #      xlab = colnames(x$df.fitted)[3], 
      #      # ylab = "vol",
      #      main = "sqrt(tau * g) and sqrt(tau) in red", sub = "",
      #      ylim = c(y_min, y_max)
      # )
      # lines(x = df_plot[, 1],
      #       y = sqrt(df_plot$tau),
      #       col = "red")
      par(mfrow = c(1, 1))
      y_max = max( max((x$g*x$tau),na.rm = TRUE) , max((x$tau),na.rm = TRUE) )
      y_min = min( min((x$g*x$tau),na.rm = TRUE) , min((x$tau),na.rm = TRUE) )
      # plot(x = df_plot[, 1], y = ((x$par[['k']])^((df_plot$e)^(-1)))*sqrt(df_plot$g*df_plot$tau),
      plot(x = x$df.fitted[,3], y = (x$g*x$tau),
           type = "l",
           xlab = colnames(x$df.fitted)[3], 
           # ylab = "vol",
           main = "(tau * g) and (tau) in red", sub = "",
           ylim = c(y_min, y_max)
      )
      lines(x = x$df.fitted[,3],
            y = (x$tau),
            col = "red")
    }
    
  }
}

#计算误差
calculate_error <- function( modwt_df, vector, method = 'mean'){
  n1 <- nrow(modwt_df)
  vector <- vector[!is.na(vector)]
  n0 <- length(vector)
  real_volatility <- c('rv','RBV','MedRV','RVKernel')
  in_sample_error=data.frame(colnames(
    c("rv_MAE", "rv_MSE", "RBV_MAE","RBV_MSE", "MedRV_MAE", "MedRV_MSE","RVKernel_MAE", "RVKernel_MSE")))
  if ( method == 'mean'){
    new_column <- c( 
      rv_MAE = mean(abs((modwt_df[[real_volatility[1]]][(n1-n0+1):n1]-(vector))))
      , rv_MSE = mean(abs((modwt_df[[real_volatility[1]]][(n1-n0+1):n1]-(vector))^2))
      , rv_HMAE = mean( abs (1- ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ) )
      , rv_HMSE = mean( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ))^2 )  
      , rv_R2log = mean( ( log( ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ) )^2 )  
      , RBV_MAE = mean(abs((modwt_df[[real_volatility[2]]][(n1-n0+1):n1]-(vector))))
      , RBV_MSE = mean(abs((modwt_df[[real_volatility[2]]][(n1-n0+1):n1]-(vector))^2))
      , RBV_HMAE = mean( abs (1- ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ) )  
      , RBV_HMSE = mean( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ))^2 )  
      , RBV_R2log = mean( ( log( ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ) )^2 ) 
      , MedRV_MAE = mean(abs((modwt_df[[real_volatility[3]]][(n1-n0+1):n1]-(vector))))
      , MedRV_MSE = mean(abs((modwt_df[[real_volatility[3]]][(n1-n0+1):n1]-(vector))^2))
      , MedRV_HMAE = mean( abs (1- ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ) )  
      , MedRV_HMSE = mean( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ))^2 )  
      , MedRV_R2log = mean( ( log( ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ) )^2 ) 
      , RVKernel_MAE = mean(abs((modwt_df[[real_volatility[4]]][(n1-n0+1):n1]-(vector))))
      , RVKernel_MSE = mean(abs((modwt_df[[real_volatility[4]]][(n1-n0+1):n1]-(vector))^2))
      , RVKernel_HMAE = mean( abs (1- ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ) )  
      , RVKernel_HMSE = mean( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ))^2 )  
      , RVKernel_R2log = mean( ( log( ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ) )^2 ) 
    )
  }

 if ( method == 'detail'){
   new_column <- list( 
     rv_MAE = (abs((modwt_df[[real_volatility[1]]][(n1-n0+1):n1]-(vector))))
     , rv_MSE = (abs((modwt_df[[real_volatility[1]]][(n1-n0+1):n1]-(vector))^2))
     , rv_HMAE = ( abs (1- ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ) )
     , rv_HMSE = ( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ))^2 )  
     , rv_R2log = ( ( log( ( vector ) / ( modwt_df[[real_volatility[1]]][(n1-n0+1):n1]) ) )^2 )  
     , RBV_MAE = (abs((modwt_df[[real_volatility[2]]][(n1-n0+1):n1]-(vector))))
     , RBV_MSE = (abs((modwt_df[[real_volatility[2]]][(n1-n0+1):n1]-(vector))^2))
     , RBV_HMAE = ( abs (1- ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ) )  
     , RBV_HMSE = ( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ))^2 )  
     , RBV_R2log = ( ( log( ( vector ) / ( modwt_df[[real_volatility[2]]][(n1-n0+1):n1]) ) )^2 ) 
     , MedRV_MAE = (abs((modwt_df[[real_volatility[3]]][(n1-n0+1):n1]-(vector))))
     , MedRV_MSE = (abs((modwt_df[[real_volatility[3]]][(n1-n0+1):n1]-(vector))^2))
     , MedRV_HMAE = ( abs (1- ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ) )  
     , MedRV_HMSE = ( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ))^2 )  
     , MedRV_R2log = ( ( log( ( vector ) / ( modwt_df[[real_volatility[3]]][(n1-n0+1):n1]) ) )^2 ) 
     , RVKernel_MAE = (abs((modwt_df[[real_volatility[4]]][(n1-n0+1):n1]-(vector))))
     , RVKernel_MSE = (abs((modwt_df[[real_volatility[4]]][(n1-n0+1):n1]-(vector))^2))
     , RVKernel_HMAE = ( abs (1- ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ) )  
     , RVKernel_HMSE = ( abs ( (1- ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ))^2 )  
     , RVKernel_R2log = ( ( log( ( vector ) / ( modwt_df[[real_volatility[4]]][(n1-n0+1):n1]) ) )^2 ) 
   )
 }
  return( as.data.frame(t(new_column)) )
}

#该函数为重新计算pot模型的似然函数值
calculate_llh <- function(df,x,y,low.freq,threshold,model, K){
  #初始化参数
  ret <- df[[x]]
  threshold_value <- threshold
  k <- c(1,1)
  par <- model[['par']]
  omega <- 1-par[['alpha']]-par[['beta']]-par[['gamma']]
  alpha <- par[['alpha']]
  beta <- par[['beta']]
  gamma <- par[['gamma']]
  m <- par[['m']]
  theta <- par[['theta']]
  Pot_p1 <- c(par[['Pot_p1_0']],par[['Pot_p1_1']])
  Pot_p2 <- c(par[['Pot_p2_0']],par[['Pot_p2_1']])
  Pot_p3 <- c(par[['Pot_p3_0']],par[['Pot_p3_1']])
  w1 <- 1
  w2 <- par[['w2']]
  tau <- model[['tau']]
  g <- model[['g']]
  
  
  ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood # is.na返回逻辑值trueorflase，1和0，which.min返回最小值
  # Pot_y <- Pot_y[which.min(is.na(tau)):length(tau)]
  g <- g[which.min(is.na(tau)):length(tau)]
  tau <- tau[which.min(is.na(tau)):length(tau)]  # 因此0的最小下标就是tau有值的第一个位置的索引
  
  sigma <- tau*g
  ret <- ret/sqrt(sigma)
  Pot_y = ret
  Pot_y[ret>0] = Pot_y[ret>0]-sort(ret)[floor(length(ret)*threshold_value)]
  Pot_y[ret<0] = Pot_y[ret<0]-sort(ret)[floor(length(ret)*(1-threshold_value))]
  Pot_y[ret>0] = ifelse(Pot_y[ret>0]>0,Pot_y[ret>0],0)
  Pot_y[ret<0] = ifelse(Pot_y[ret<0]<0,abs(Pot_y[ret<0]),0)
  threshold_value = c(sort(ret)[floor(length(ret)*threshold_value)], -sort(ret)[floor(length(ret)*(1-threshold_value))])
  yt = sign(Pot_y)
  sigma1<-rep(1,length(ret))
  
  if(length(threshold_value)>1){
    e_0<-calculate_shape_value(b1 = Pot_p1[1],b2 = Pot_p2[1],b3 = Pot_p3[1],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
    e_1<-calculate_shape_value(b1 = Pot_p1[2],b2 = Pot_p2[2],b3 = Pot_p3[2],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
    
    penalty0 <- ifelse( ( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) )>0
                        , log( 1-k[1]*(threshold_value[1]/sigma1)^(-e_0)-k[2]*(threshold_value[2]/sigma1)^(-e_1) ),-1e6 )
    penalty1 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                        , ( log(k[1])+log(e_0)+e_0*log(sigma1)-(e_0+1)*log(threshold_value[1]+Pot_y) )
                        , -1e6 )
    penalty2 <- ifelse( ( ( k[1]*( (threshold_value[1]) /sigma1)^(-e_0) ) + k[2]*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                        , ( log(k[2])+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]+Pot_y) )
                        , -1e6 )
    return(
      -sum(-(((yt-1)/-1)*penalty0
             +yt*(ret>0)*(penalty1)
             +yt*(ret<0)*(penalty2)))
    )#注意区分小k和大K
  }
  
}

# 超多频段外生变量滚动预测
rolling_predict_volatility <- function(i){
  #滚动预测
  predict_volatility <- numeric( length(order_set) )
  tau_forecast <- numeric( length(order_set) )
  uptail_parameters <- numeric( length(order_set) )
  downtail_parameters <- numeric( length(order_set) )
  tail_parameters <- numeric( length(order_set) )
  
  DataFrame=df[(i-window_length+1):i,]#一定要括号起来,不然是错误的
  for (j in 1:length(order_set)){
    numbers <- str_extract_all(order_set[j], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
    if(grepl('two_side_Pot',order_set[j])){
      fit_model = Pot_fit_garch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                , pot_k = TRUE, two_side = TRUE, gamma = TRUE )
      e_0 <- calculate_shape_value( fit_model$par[['Pot_p1_0']], fit_model$par[['Pot_p2_0']], fit_model$par[['Pot_p3_0']]
                                    , c(fit_model$df.fitted$return,0) )# 多一个长度是为了预测t+1时刻的尾部参数
      e_1 <- calculate_shape_value( fit_model$par[['Pot_p1_1']], fit_model$par[['Pot_p2_1']], fit_model$par[['Pot_p3_1']]
                                    , c(fit_model$df.fitted$return,0) )
      predict_volatility[j] <- predict(fit_model, c(1) ,method = method)
      tau_forecast[j] <- fit_model$tau.forecast
      uptail_parameters[j] <- tail(e_0,1)
      downtail_parameters[j] <- tail(e_1,1)
      

    }else if(grepl('one_side_Pot',order_set[j])){
      fit_model = Pot_fit_garch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                , pot_k = TRUE, two_side = FALSE, gamma = TRUE )
      e <- calculate_shape_value( fit_model$par[['Pot_p1']], fit_model$par[['Pot_p2']], fit_model$par[['Pot_p3']]
                                  , c(fit_model$df.fitted$return,0) )# 多一个长度是为了预测t+1时刻的尾部参数
      predict_volatility[j] <- predict(fit_model, c(1) ,method = method)
      tau_forecast[j] <- fit_model$tau.forecast
      tail_parameters[j] <- tail(e,1)
      
    }else{
      fit_model=fit_mfgarch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[1]))
                            , low.freq = low.freq, K = K,distribution = str_extract(order_set[j], "^[a-z]+" ))
      predict_volatility[j] <- predict(fit_model, c(1) ,method = method)      }
      tau_forecast[j] <- fit_model$tau.forecast
    
  }
  return( list( predict_volatility = predict_volatility , tau_forecast = tau_forecast , uptail_parameters = uptail_parameters
                , downtail_parameters = downtail_parameters, tail_parameters = tail_parameters) ) 
}

# 计算超阈值部分的累计分布
{
Pot_y_cdf <- function( ret, parameters, threshold){
  if( length(threshold) == 2 ){
    threshold_0 <- threshold[["threshold_0"]]
    threshold_1 <- threshold[["threshold_1"]]
    #计算阈值
    # threshold_0 <- sort(ret[ret>0])[floor(length(ret[ret>0])*threshold_0)]
    # threshold_1 <- abs(sort(ret[ret<0])[floor(length(ret[ret<0])*(1-threshold_1))])
    threshold_0 <- sort(ret)[floor(length(ret)*threshold_0)]
    threshold_1 <- abs(sort(ret)[floor(length(ret)*(1-threshold_1))]) # 由于threshold_1的值在代码中要被用于计算分布，因此采用绝对值方便书写，避免-+号的混乱
    threshold <- c( threshold_0 = threshold_0, threshold_1  = threshold_1) 
    # 计算Pot_y
    Pot_y = ret
    Pot_y[ret>0] = Pot_y[ret>0]-threshold_0
    Pot_y[ret<0] = Pot_y[ret<0]+threshold_1
    Pot_y[ret>0] = ifelse( Pot_y[ret>0]>0, Pot_y[ret>0], 0)
    Pot_y[ret<0] = ifelse( Pot_y[ret<0]<0, Pot_y[ret<0], 0)
    # 计算概率
    return(Pot_distribution( Pot_y, parameters,threshold))
  }else if( length(threshold) == 1 ){
    # 计算阈值
    threshold <- sort(ret)[floor(length(ret)*threshold)]
    # 计算Pot_y
    Pot_y <- (ret-threshold)
    Pot_y = ifelse(Pot_y>0,Pot_y,0)
    #计算阈值
    return(Pot_distribution( Pot_y, parameters,threshold))
  }
  
}

Pot_distribution <- function( Pot_y, parameters, threshold){
  if( length(threshold) ==2){
    #初始化参数
    sigma <- parameters[["sigma"]]
    Pot_y <- Pot_y[ which.max( !is.na(sigma) ): length(sigma) ]
    e_0 <- parameters[["e_0"]][ which.max( !is.na(sigma) ): length(sigma) ]
    e_1 <- parameters[["e_1"]][ which.max( !is.na(sigma) ): length(sigma) ]
    k_0 <- parameters[["k_0"]]
    k_1 <- parameters[["k_1"]]
    sigma <- sigma[!is.na(sigma)]
    sigma <-sqrt(sigma)
    threshold_0 <- threshold[["threshold_0"]]
    threshold_1 <- threshold[["threshold_1"]]
    # 计算概率
    Pot_yy <- Pot_y
    Pot_yy[Pot_y < 0] <- k_1*( ( -Pot_yy[Pot_y < 0]+threshold_1 )/sigma[Pot_y < 0])^(-e_1[Pot_y < 0])
    Pot_yy[Pot_y == 0] <- 1-k_0*( ( threshold_0 )/sigma[Pot_y == 0])^(-e_0[Pot_y == 0])
    Pot_yy[Pot_y > 0] <- 1-k_0*( ( Pot_yy[Pot_y > 0]+threshold_0 )/sigma[Pot_y > 0])^(-e_0[Pot_y > 0])
    return(Pot_yy)
  }else if( length(threshold) == 1){
    
    #初始化参数
    sigma <- parameters[["sigma"]]
    Pot_y <- Pot_y[ which.max( !is.na(sigma) ): length(sigma) ]
    e <- parameters[["e"]][ which.max( !is.na(sigma) ): length(sigma) ]
    k <- parameters[["k"]]
    sigma <- sigma[!is.na(sigma)]
    sigma <-sqrt(sigma)
    threshold <- threshold
    # 计算概率
    Pot_yy <- Pot_y
    Pot_yy[Pot_y == 0] <- 1-k*( ( threshold )/sigma[Pot_y == 0] )^( -e[Pot_y == 0] )
    Pot_yy[Pot_y > 0] <- 1-k*( ( Pot_yy[Pot_y > 0]+threshold )/sigma[Pot_y > 0])^(-e[Pot_y > 0])
    return(Pot_yy)
  }

}
}
# 求解双边pot_cdf的逆分布、或是收益率
two_side_reverse_Pot_distribution <- function( u, parameters, ret, threshold){
  # 计算阈值
  threshold_0 <- threshold[["threshold_0"]]
  threshold_1 <- threshold[["threshold_1"]]
  threshold_0 <- sort(ret)[floor(length(ret)*threshold_0)]
  threshold_1 <- abs(sort(ret)[floor(length(ret)*(1-threshold_1))]) # 由于threshold_1的值在代码中要被用于计算分布，因此采用绝对值方便书写，避免-+号的混乱
  threshold <- c( threshold_0 = threshold_0, threshold_1  = threshold_1) 
  ff <- function(theta) { ( 100*( Pot_distribution( theta, parameters, threshold)-u))^2 }
  #乘以10是为了扩大函数，这样子找最小值会更加精确
  result <-  optim(0.5, function(x) ff(x))$par
  return(result)
}
# 求解单边pot_cdf的逆分布、或是收益率
one_side_reverse_Pot_distribution <- function( u, parameters, ret, threshold){
  # 计算阈值
  threshold <- threshold
  threshold <- sort(ret)[floor(length(ret)*threshold)]
  ff <- function(theta) { ( 100*( Pot_distribution( theta, parameters, threshold)-u))^2 }
  #乘以10是为了扩大函数，这样子找最小值会更加精确
  result <-  optim(0.5, function(x) ff(x))$par
  return(result)
}

# 计算的残差对应的累积分布
residuals_distribution <- function(x,parameters,data,threshold){
  #初始化参数
  sigma <- parameters[["sigma"]]
  e_0 <- parameters[["e_0"]]
  e_1 <- parameters[["e_1"]]
  k_0 <- parameters[["k_0"]]
  k_1 <- parameters[["k_1"]]
  sigma <-sqrt(sigma)
  threshold_0 <- threshold[["threshold_0"]]
  threshold_1 <- threshold[["threshold_1"]]
  #计算阈值
  threshold_0 <- sort(data)[floor(length(data)*threshold_0)]
  threshold_1 <- abs(sort(data)[floor(length(data)*(1-threshold_1))]) # 由于threshold_1的值在代码中要被用于计算分布，因此采用绝对值方便书写，避免-+号的混乱
  #定义截断的经验分布函数，也就是小于阈值下的经验分布函数
  empirical_distribution <- function(x,data,threshold_0,threshold_1){
    return( sum((data<x)&(data<threshold_0)&(data>(-threshold_1)) )/sum( (data<threshold_0)&(data>(-threshold_1)) ) )
  }
  #
  
  
  # if(x>((1/k)^(-1/e)*threshold))
  if(x>threshold_0)
  {
    # print((1-k*(threshold/sigma)^(-e))*empirical_distribution(x,data,threshold))
    return(
      (1-k_0*(threshold_0/sigma)^(-e_0)-k_1*(threshold_1/sigma)^(-e_1))
      *empirical_distribution(x,data,threshold_0,threshold_1)
      # +k_0*( (threshold_0/sigma)^(-e_0) )*( 1-k_0*(x/threshold_0)^(-e_0) )
      +k_0*( (threshold_0/sigma)^(-e_0) )*( 1-(x/threshold_0)^(-e_0) )
      +k_1*( (threshold_1/sigma)^(-e_1) )
    )
  }else if( (x<threshold_0)&((-threshold_1)<x) ){
    return(
      (1-k_0*(threshold_0/sigma)^(-e_0)-k_1*(threshold_1/sigma)^(-e_1) )
      *empirical_distribution(x,data,threshold_0,threshold_1)
      +k_1*( (threshold_1/sigma)^(-e_1) )
    )
  }else if(x<(-threshold_1)){
    # return( k_1*( (threshold_1/sigma)^(-e_1) )*( 1-( 1-k_1*(-x/threshold_1)^(-e_1) ) ) )
    return( k_1*( (threshold_1/sigma)^(-e_1) )*( 1-( 1-(-x/threshold_1)^(-e_1) ) ) )
    
  }
  
}

rolling_predict_garch <- function(i) {
  # 截取窗口数据
  DataFrame <- df[(i - window_length + 1):i,]
  # 创建一个预测集合
  predict_set <- list()
  for (j in 1:length(order_set)) {
    predict_set[[order_set[j]]] <- as.numeric(NA)
  }
  # 拟合不同模型
  for (j in 1:length(order_set)) {
    numbers <- str_extract_all(order_set[j], "\\d+\\.?\\d*")[[1]]
    if (grepl('Pot', order_set[j])) {
      fit_model <- Pot_fit_garch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                 low.freq = low.freq, K = K, threshold_value = as.numeric(numbers[1]),
                                 pot_k = FALSE, two_side = TRUE, gamma = TRUE) 
      predict_set[[order_set[j]]] <- predict(fit_model, c(1), method = 'normal')
    } else {
      fit_model <- fit_mfgarch(data = DataFrame, y = "return", x = paste0("rv", as.numeric(numbers[1])),
                               low.freq = low.freq, K = K, distribution = str_extract(order_set[j], "^[a-z]+"))
      predict_set[[order_set[j]]] <- predict(fit_model, c(1), method = 'normal')
    }
  }
  
  return( predict_set )
}
