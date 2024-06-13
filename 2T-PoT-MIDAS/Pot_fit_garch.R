#' This function estimates a multiplicative mixed-frequency GARCH model. For the sake of numerical stability, it is best to multiply log returns by 100.
#' @param data data frame containing a column named date of type 'Date'.
#' @param y name of high frequency dependent variable in df.
#' @param x covariate employed in mfGARCH.
#' @param K an integer specifying lag length K in the long-term component.
#' @param low.freq a string of the low frequency variable in the df.
#' @param var.ratio.freq specify a frequency column on which the variance ratio should be calculated.
#' @param gamma if TRUE, an asymmetric GJR-GARCH is used as the short-term component. If FALSE, a simple GARCH(1,1) is employed.
#' @param weighting specifies the weighting scheme employed in the long-term component. Options are "beta.restricted" (default) or "beta.unrestricted"
#' @param x.two optional second covariate
#' @param K.two lag lgenth of optional second covariate
#' @param low.freq.two low frequency of optional second covariate
#' @param weighting.two specifies the weighting scheme employed in the optional second long-term component. Currently, the only option is "beta.restricted"
#' @param multi.start if TRUE, optimization is carried out with multiple starting values
#' @param control a list
#' @keywords fit_mfgarch
#' @export
#'
#' @return A list of class mfGARCH with letters and numbers.
#' \itemize{
#' \item par - vector of estimated parameters
#' \item rob.std.err - sandwich/HAC-type standard errors
#' \item broom.mgarch - a broom-like data.frame with entries
#' 1) estimate: column of estimated parameters
#' 2) rob.std.err - sandwich/HAC-type standard errors
#' 3) p.value - p-values derived from sandwich/HAC-type standard errors
#' 4) opg.std.err - Bollerslev-Wooldrige/OPG standard errors for GARCH processes
#' 5) opg.p.value - corresponding alternative p-values
#' \item tau - fitted long-term component
#' \item g - fitted short-term component
#' \item df.fitted - data frame with fitted values and residuals
#' \item K - chosen lag-length in the long-term component
#' \item weighting.scheme - chosen weighting scheme
#' \item llh - log-likelihood value at estimated parameter vector
#' \item bic - corresponding BIC value
#' \item y - dependent variable y
#' \item optim - output of the optimization routine
#' \item K.two - lag-lenth of x.two if two covariates are employed
#' \item weighting.scheme.two - chosen weighting scheme of x.two (if K.two != NULL)
#' \item tau.forecast - one-step ahead forecast of the long-term component
#' \item variance.ratio - calculated variance ratio
#' \item est.weighting - estimated weighting scheme
#' \item est.weighting.two - estimated weighting scheme of x.two (if K.two != NULL)
#' }
#'
#' @importFrom numDeriv jacobian
#' @importFrom stats nlminb
#' @importFrom numDeriv hessian
#' @importFrom stats constrOptim
#' @importFrom stats na.exclude
#' @importFrom stats optim
#' @importFrom stats pnorm
#' @importFrom stats var
#' @importFrom stats aggregate
#' @importFrom numDeriv jacobian
#' @importFrom maxLik maxLik
#' @importFrom utils tail
#' @examples
#' \dontrun{
#' fit_mfgarch(data = df_financial, y = "return", x = "nfci", low.freq = "week", K = 52)
#' fit_mfgarch(data = df_mfgarch, y = "return", x = "nfci", low.freq = "year_week", K = 52,
#' x.two = "dindpro", K.two = 12, low.freq.two = "year_month", weighting.two = "beta.restricted")
#' }

Pot_fit_garch <- function(data, y, x = NULL, K = NULL, low.freq , var.ratio.freq = NULL, gamma = FALSE, weighting = "beta.restricted"
                          , x.two = NULL, K.two = NULL, low.freq.two = NULL, weighting.two = NULL, multi.start = FALSE
                          , control = list(par.start = NULL), threshold_value = 0.9, methods='Pot', pot_k = TRUE, two_side = FALSE) {
  
  print("For ensuring numerical stability of the parameter optimization and inversion of the Hessian, it is best to multiply log returns by 100.")

  data <- data[order(data$date), ]
  date_backup <- data[["date"]]

  if (is.null(var.ratio.freq)) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }
  
  low_freq_backup <- data[, low.freq]
  if (x != "date") {
    if (is.null(x.two)) {
      df_llh <- data[, c(y, x, low.freq)]
      df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
    } else {
      low_freq.two_backup <- data[, low.freq.two]
      if (low.freq != low.freq.two) { # if they are different, both have to be included in df_llh
        df_llh <- data[, c(y, x, low.freq, x.two, low.freq.two)]
        df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
        df_llh[, low.freq.two] <- as.integer(unlist(df_llh[ , low.freq.two]))
      } else { # else, the low.freq column is needed only once
        df_llh <- data[, c(y, x, low.freq, x.two)]
        df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
      }
    }
  }
  
  g_zero <- var(unlist(data[[y]]))
  ret <- data[[y]]
  # ret <- ret-mean(ret)去均值效果不好，POT研究的是收益尾部，原始收益不要进行处理
  # Parameter estimation
  # if K > 0 we get the covariate series
  covariate <- unlist(unique(data[c(low.freq, x)])[,x])
  #提取出低频的数据
  
  # 计算出阈值
  if (K > 1) {
    if (methods=='Pot') {
      if(two_side==TRUE){
        Pot_y=ret
        Pot_y[ret>0] = Pot_y[ret>0]-sort(ret)[floor(length(ret)*threshold_value)]
        Pot_y[ret<0] = Pot_y[ret<0]-sort(ret)[floor(length(ret)*(1-threshold_value))]
        Pot_y[ret>0] = ifelse(Pot_y[ret>0]>0,Pot_y[ret>0],0)
        Pot_y[ret<0] = ifelse(Pot_y[ret<0]<0,abs(Pot_y[ret<0]),0)
        if ( is.null(K.two) == TRUE) {
          if(pot_k==FALSE){
            if(gamma==TRUE){
              lf <- function(p) {
                llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                           , k = c(1,1)
                           , omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2
                           , alpha = p["alpha"]
                           , beta = p["beta"], gamma = p["gamma"]
                           , m = p["m"], theta = p["theta"]
                           , w1 = 1, w2 = p["w2"]
                           # , threshold_value=c(sort(ret[ret>0])[floor(length(ret[ret>0])*threshold_value)]
                           #                     ,abs(sort(ret[ret<0])[floor(length(ret[ret<0])*(1-threshold_value))]))
                           , threshold_value=c(sort(ret)[floor(length(ret)*threshold_value)]
                                               ,abs(sort(ret)[floor(length(ret)*(1-threshold_value))]))
                           , Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"])
                           , Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"])
                           , Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"])
                           , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                           , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
              }
              
              par.start <- c(#mu = 0,
                alpha = 0.0005, beta = 0.91, gamma = 0.0005,
                m = 0, theta = -0.01, w2 = 3, 
                Pot_p1_0 = 0.1, Pot_p2_0 = 0.6, Pot_p3_0 = 0.1,
                Pot_p1_1 = 0.1, Pot_p2_1 = 0.6, Pot_p3_1 = 0.1)
              
              # ui.opt <- rbind(c(-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0),
              #                 c(0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
              #                 c(0,  0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
              #                 c(1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
              #                 c(0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0,-1,-1, 0, 0, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
              #                 c(0,  0,  0, 0, 0, 0, 0, 0, 0, 0,-1,-1))
              # 
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0,-1, 0,-0.99999999,-1, 0,-0.99999999)
              ui.opt <- rbind(c(-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0),
                              c(0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                              c(0,  0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                              c(1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                              c(0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                              c(0,  0,  0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                              c(0,  0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                              c(0,  0,  0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                              c(0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

              ci.opt <- c(-0.99999999, 0, 1, 0, 0,-1, 0,-1, 0)
            }else if(gamma==FALSE){
              lf <- function(p) {
                llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                           , k = c(1,1)
                           , omega = 1 - p["alpha"] - p["beta"]# - p["gamma"]/2
                           , alpha = p["alpha"]
                           , beta = p["beta"], gamma = 0
                           , m = p["m"], theta = p["theta"]
                           , w1 = 1, w2 = p["w2"]
                           , threshold_value=c(sort(ret)[floor(length(ret)*threshold_value)]
                                               ,abs(sort(ret)[floor(length(ret)*(1-threshold_value))]))
                           , Pot_p1 = c(p["Pot_p1"],p["Pot_p1"])
                           , Pot_p2 = c(p["Pot_p2"],p["Pot_p2"])
                           , Pot_p3 = c(p["Pot_p3"],p["Pot_p3"])
                           , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                           , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
              }
              
              par.start <- c(#mu = 0,
                alpha = 0.0005, beta = 0.91,# gamma = 0.0005,
                m = 0, theta = -0.01, w2 = 3, 
                Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
              
              # ui.opt <- rbind(c(-1,-1, 0, 0, 0, 0, 0, 0),
              #                 c(0,  1, 0, 0, 0, 0, 0, 0),
              #                 c(0,  0, 0, 0, 1, 0, 0, 0),
              #                 c(1,  0, 0, 0, 0, 0, 0, 0),
              #                 c(0,  1, 0, 0, 0, 0, 0, 0),
              #                 c(0,  0, 0, 0, 0, 0,-1, 0),
              #                 c(0,  0, 0, 0, 0, 0, 1, 0),
              #                 c(0,  0, 0, 0, 0, 0,-1,-1))
              # 
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0,-1, 0,-0.99999999)
              ui.opt <- rbind(c(-1,-1, 0, 0, 0, 0, 0, 0),
                              c(0,  1, 0, 0, 0, 0, 0, 0),
                              c(0,  0, 0, 0, 1, 0, 0, 0),
                              c(1,  0, 0, 0, 0, 0, 0, 0),
                              c(0,  1, 0, 0, 0, 0, 0, 0),
                              c(0,  0, 0, 0, 0, 0,-1, 0),
                              c(0,  0, 0, 0, 0, 0, 1, 0))
              
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0,-1, 0)
              ci.opt <- c(-0.99999999,0, 1, 0, 0,-1, 0)
            }
          }
          else if(pot_k==TRUE){
            if(gamma==TRUE){
              if(weighting == "beta.unrestricted"){
                lf <- function(p) {
                  
                  llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                             , k = c(p["k_0"], p["k_1"])
                             , omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2
                             , alpha = p["alpha"]
                             , beta = p["beta"], gamma = p["gamma"]
                             , m = p["m"], theta = p["theta"]
                             , w1 = p["w1"], w2 = p["w2"]
                             , threshold_value=c(sort(ret)[floor(length(ret)*threshold_value)]
                                                 ,abs(sort(ret)[floor(length(ret)*(1-threshold_value))]))
                             , Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"])
                             , Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"])
                             , Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"])
                             , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                             , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
                }
                
                par.start <- c(#mu = 0,
                  k_0= 0.01,k_1=0.01,
                  alpha = 0.0005, beta = 0.91, gamma = 0.0005,
                  m = 0, theta = -0.01, w1 = 1.0001,w2 = 3,
                  Pot_p1_0 = 0.1, Pot_p2_0 = 0.6, Pot_p3_0 = 0.1,
                  Pot_p1_1 = 0.1, Pot_p2_1 = 0.6, Pot_p3_1 = 0.1)
                # ui.opt <- rbind(c(0, 0,-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0,-1,-1, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1))
                # 
                # ci.opt <- c(-0.99999999, 0.8, 1, 1, 0, 0, -1, 0, -0.99999999, -1, 0, -0.99999999)
                ui.opt <- rbind(c(0, 0,-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
                
                # ci.opt <- c(-0.99999999, 0.8, 1, 1, 0, 0, -1, 0, -1, 0)
                ci.opt <- c(-0.99999999, 0, 1, 1, 0, 0, -1, 0, -1, 0)
              }else{
                lf <- function(p) {
                  
                  llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                             , k = c(p["k_0"], p["k_1"])
                             , omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2
                             , alpha = p["alpha"]
                             , beta = p["beta"], gamma = p["gamma"]
                             , m = p["m"], theta = p["theta"]
                             , w1 = 1, w2 = p["w2"]
                             , threshold_value=c(sort(ret)[floor(length(ret)*threshold_value)]
                                                 ,abs(sort(ret)[floor(length(ret)*(1-threshold_value))]))
                             , Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"])
                             , Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"])
                             , Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"])
                             , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                             , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
                }
                
                par.start <- c(#mu = 0,
                  k_0= 0.01,k_1=0.01,
                  alpha = 0.0005, beta = 0.91, gamma = 0.0005,
                  m = 0, theta = -0.01, w2 = 3,
                  Pot_p1_0 = 0.1, Pot_p2_0 = 0.6, Pot_p3_0 = 0.1,
                  Pot_p1_1 = 0.1, Pot_p2_1 = 0.6, Pot_p3_1 = 0.1)
                # ui.opt <- rbind(c(0, 0,-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0,-1,-1, 0, 0, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                #                 c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,-1,-1))
                # 
                # ci.opt <- c(-0.99999999,0.8, 1, 0, 0,-1, 0,-0.99999999,-1, 0,-0.99999999)
                ui.opt <- rbind(c(0, 0,-1,-1,-1/2,0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                                c(0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
                
                # ci.opt <- c(-0.99999999,0.8, 1, 0, 0,-1, 0, -1, 0)
                ci.opt <- c(-0.99999999,0, 1, 0, 0,-1, 0, -1, 0)
                
              }
          }else if(gamma==FALSE){
              lf <- function(p) {
                
                llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                           , k = c(p["k_0"], p["k_1"])
                           , omega = 1 - p["alpha"] - p["beta"]# - p["gamma"]/2
                           , alpha = p["alpha"]
                           , beta = p["beta"], gamma = 0
                           , m = p["m"], theta = p["theta"]
                           , w1 = 1, w2 = p["w2"]
                           , threshold_value=c(sort(ret)[floor(length(ret)*threshold_value)]
                                               ,abs(sort(ret)[floor(length(ret)*threshold_value)]))
                           , Pot_p1 = c(p["Pot_p1"],p["Pot_p1"])
                           , Pot_p2 = c(p["Pot_p2"],p["Pot_p2"])
                           , Pot_p3 = c(p["Pot_p3"],p["Pot_p3"])
                           , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                           , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
              }
              
              par.start <- c(#mu = 0,
                k_0= 0.01,k_1=0.01,
                alpha = 0.0005, beta = 0.91,# gamma = 0.0005,
                m = 0, theta = -0.01, w2 = 3, Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
              # ui.opt <- rbind(c(0, 0,-1, -1,  0, 0, 0, 0, 0, 0),
              #                 c(0, 0, 0,  1,  0, 0, 0, 0, 0, 0),
              #                 c(0, 0, 0,  0,  0, 0, 1, 0, 0, 0),
              #                 c(0, 0, 1,  0,  0, 0, 0, 0, 0, 0),
              #                 c(0, 0, 0,  1,  0, 0, 0, 0, 0, 0),
              #                 c(0, 0, 0,  0,  0, 0, 0, 0,-1, 0),
              #                 c(0, 0, 0,  0,  0, 0, 0, 0, 1, 0),
              #                 c(0, 0, 0,  0,  0, 0, 0, 0, -1, -1))
              # 
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0,-0.99999999)
              ui.opt <- rbind(c(0, 0,-1, -1,  0, 0, 0, 0, 0, 0),
                              c(0, 0, 0,  1,  0, 0, 0, 0, 0, 0),
                              c(0, 0, 0,  0,  0, 0, 1, 0, 0, 0),
                              c(0, 0, 1,  0,  0, 0, 0, 0, 0, 0),
                              c(0, 0, 0,  1,  0, 0, 0, 0, 0, 0),
                              c(0, 0, 0,  0,  0, 0, 0, 0,-1, 0),
                              c(0, 0, 0,  0,  0, 0, 0, 0, 1, 0))
              
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0)
              ci.opt <- c(-0.99999999,0, 1, 0, 0, -1, 0)
            }
 
          }
        }
      }
      else if(two_side==FALSE){
        # Pot_y = ret-sort(ret)[floor(length(ret)*threshold_value)]
        Pot_y = (ret-sort(ret)[floor(length(ret)*threshold_value)])
        Pot_y = ifelse(Pot_y>0,Pot_y,0)
        if (is.null(K.two) == TRUE) {
          if(pot_k==FALSE){
             if( gamma == TRUE ){
              lf <- function(p) {
                llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                           , k = 1
                           # , omega = 1 - p["alpha"] - p["beta"]# - p["gamma"]/2
                           , omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2
                           , alpha = p["alpha"]
                           , beta = p["beta"], gamma = p["gamma"]
                           , m = p["m"], theta = p["theta"]
                           , w1 = 1, w2 = p["w2"]
                           , threshold_value=sort(ret)[floor(length(ret)*threshold_value)]
                           , Pot_p1 = p["Pot_p1"], Pot_p2 = p["Pot_p2"], Pot_p3 = p["Pot_p3"]
                           , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                           , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
              }
              
              par.start <- c(# mu = 0,
                alpha = 0.0005, beta = 0.91, gamma = 0.0005,
                m = 0, theta = -0.01, w2 = 3, Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
              
              # ui.opt <- rbind(c(-1, -1, 0, 0, 0, 0, 0, 0, 0),
              #                 c( 0,  1, 0, 0, 0, 0, 0, 0, 0),
              #                 c( 0,  0, 0, 0, 0, 1, 0, 0, 0),
              #                 c( 1,  0, 0, 0, 0, 0, 0, 0, 0),
              #                 c( 0,  1, 0, 0, 0, 0, 0, 0, 0),
              #                 c( 0,  0, 0, 0, 0, 0, 0,-1, 0),
              #                 c( 0,  0, 0, 0, 0, 0, 0, 1, 0),
              #                 c( 0,  0, 0, 0, 0, 0, 0, -1, -1))
              # 
              # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0,-0.99999999)
              ui.opt <- rbind(c(-1, -1, 0, 0, 0, 0, 0, 0, 0),
                              c( 0,  1, 0, 0, 0, 0, 0, 0, 0),
                              c( 0,  0, 0, 0, 0, 1, 0, 0, 0),
                              c( 1,  0, 0, 0, 0, 0, 0, 0, 0),
                              c( 0,  1, 0, 0, 0, 0, 0, 0, 0),
                              c( 0,  0, 0, 0, 0, 0, 0,-1, 0),
                              c( 0,  0, 0, 0, 0, 0, 0, 1, 0))
              
              ci.opt <- c(-0.99999999, 0, 1, 0, 0, -1, 0)
              
            }else if( gamma == FALSE ){
            lf <- function(p) {
              llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                         , k = 1
                         , omega = 1 - p["alpha"] - p["beta"]# - p["gamma"]/2
                         , alpha = p["alpha"]
                         , beta = p["beta"], gamma = 0
                         , m = p["m"], theta = p["theta"]
                         , w1 = 1, w2 = p["w2"]
                         , threshold_value=sort(ret)[floor(length(ret)*threshold_value)]
                         , Pot_p1 = p["Pot_p1"], Pot_p2 = p["Pot_p2"], Pot_p3 = p["Pot_p3"]
                         , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                         , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
            }
            
            par.start <- c(# mu = 0,
              alpha = 0.0005, beta = 0.91,# gamma = 0.0005,
              m = 0, theta = -0.01, w2 = 3, Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
            
            # ui.opt <- rbind(c(-1, -1,  0, 0, 0, 0, 0, 0),
            #                 c( 0,  1,  0, 0, 0, 0, 0, 0),
            #                 c( 0,  0,  0, 0, 1, 0, 0, 0),
            #                 c( 1,  0,  0, 0, 0, 0, 0, 0),
            #                 c( 0,  1,  0, 0, 0, 0, 0, 0),
            #                 c( 0,  0,  0, 0, 0, 0,-1, 0),
            #                 c( 0,  0,  0, 0, 0, 0, 1, 0),
            #                 c( 0,  0,  0, 0, 0, 0, -1, -1))
            # 
            # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0,-0.99999999)
            ui.opt <- rbind(c(-1, -1,  0, 0, 0, 0, 0, 0),
                            c( 0,  1,  0, 0, 0, 0, 0, 0),
                            c( 0,  0,  0, 0, 1, 0, 0, 0),
                            c( 1,  0,  0, 0, 0, 0, 0, 0),
                            c( 0,  1,  0, 0, 0, 0, 0, 0),
                            c( 0,  0,  0, 0, 0, 0,-1, 0),
                            c( 0,  0,  0, 0, 0, 0, 1, 0))
            
            # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0)
              ci.opt <- c(-0.99999999,0, 1, 0, 0, -1, 0)
            }
          }
          else if(pot_k==TRUE){
              if(gamma==TRUE){
                
                lf <- function(p) {
                llh_mf_Pot(df = df_llh, y = ret,Pot_y = Pot_y
                           , x = covariate, low.freq = low.freq, mu = 0
                           , k = p["k"]
                           , omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2
                           , alpha = p["alpha"]
                           , beta = p["beta"], gamma = p["gamma"]
                           , m = p["m"], theta = p["theta"]
                           , w1 = 1, w2 = p["w2"]
                           , threshold_value=sort(ret)[floor(length(ret)*threshold_value)]
                           , Pot_p1 = p["Pot_p1"], Pot_p2 = p["Pot_p2"], Pot_p3 = p["Pot_p3"]
                           , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                           , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
              }
              
              par.start <- c(#mu = 0, 
                             k = 0.01, alpha = 0.0005, beta = 0.91, gamma = 0.0005,
                m = 0, theta = -0.01, w2 = 3, Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
              # ui.opt <- rbind(c( 0,-1, -1,-1/2, 0, 0, 0, 0, 0, 0),
              #                 c( 0, 0,  1,  0, 0, 0, 0, 0, 0, 0),
              #                 c( 0, 0,  0,  0, 0, 0, 1, 0, 0, 0),
              #                 c( 0, 1,  0,  0, 0, 0, 0, 0, 0, 0),
              #                 c( 0, 0,  1,  0, 0, 0, 0, 0, 0, 0),
              #                 c( 0, 0,  0,  0, 0, 0, 0, 0,-1, 0),
              #                 c( 0, 0,  0,  0, 0, 0, 0, 0, 1, 0),
              #                 c( 0, 0,  0,  0, 0, 0, 0, 0, -1, -1))
              ui.opt <- rbind(c( 0,-1, -1,-1/2, 0, 0, 0, 0, 0, 0),
                              c( 0, 0,  1,  0, 0, 0, 0, 0, 0, 0),
                              c( 0, 0,  0,  0, 0, 0, 1, 0, 0, 0),
                              c( 0, 1,  0,  0, 0, 0, 0, 0, 0, 0),
                              c( 0, 0,  1,  0, 0, 0, 0, 0, 0, 0),
                              c( 0, 0,  0,  0, 0, 0, 0, 0,-1, 0),
                              c( 0, 0,  0,  0, 0, 0, 0, 0, 1, 0))
            
              }
            else if(gamma==FALSE){
              lf <- function(p) {
              llh_mf_Pot(df = df_llh, y = ret,Pot_y, x = covariate, low.freq = low.freq#, mu = p["mu"]
                         , k = p["k"]
                         , omega = 1 - p["alpha"] - p["beta"]# - p["gamma"]/2
                         , alpha = p["alpha"]
                         , beta = p["beta"], gamma = 0
                         , m = p["m"], theta = p["theta"]
                         , w1 = 1, w2 = p["w2"]
                         , threshold_value=sort(ret)[floor(length(ret)*threshold_value)]
                         , Pot_p1 = p["Pot_p1"], Pot_p2 = p["Pot_p2"], Pot_p3 = p["Pot_p3"]
                         , g_zero = g_zero, K = K,x.two = NULL, K.two = NULL, theta.two = NULL
                         , low.freq.two = NULL, w1.two = NULL, w2.two = NULL)
            }
            
            par.start <- c(#mu = 0,
              k = 0.01,
              alpha = 0.0005, beta = 0.91,# gamma = 0.0005,
              m = 0, theta = -0.01, w2 = 3, Pot_p1 = 0.1, Pot_p2 = 0.6, Pot_p3 = 0.1)
            # ui.opt <- rbind(c(0,-1, -1,  0, 0, 0, 0, 0, 0),
            #                 c(0, 0,  1,  0, 0, 0, 0, 0, 0),
            #                 c(0, 0,  0,  0, 0, 1, 0, 0, 0),
            #                 c(0, 1,  0,  0, 0, 0, 0, 0, 0),
            #                 c(0, 0,  1,  0, 0, 0, 0, 0, 0),
            #                 c(0, 0,  0,  0, 0, 0, 0,-1, 0),
            #                 c(0, 0,  0,  0, 0, 0, 0, 1, 0),
            #                 c(0, 0,  0,  0, 0, 0, 0, -1, -1))
            ui.opt <- rbind(c(0,-1, -1,  0, 0, 0, 0, 0, 0),
                            c(0, 0,  1,  0, 0, 0, 0, 0, 0),
                            c(0, 0,  0,  0, 0, 1, 0, 0, 0),
                            c(0, 1,  0,  0, 0, 0, 0, 0, 0),
                            c(0, 0,  1,  0, 0, 0, 0, 0, 0),
                            c(0, 0,  0,  0, 0, 0, 0,-1, 0),
                            c(0, 0,  0,  0, 0, 0, 0, 1, 0))
            }
            
            # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0,-0.99999999)
            # ci.opt <- c(-0.99999999,0.8, 1, 0, 0, -1, 0)
            ci.opt <- c(-0.99999999,0, 1, 0, 0, -1, 0)

          }
        }
        
      }
    }
    p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)

    p.e.nlminb$value <- -p.e.nlminb$value
    
    par <- p.e.nlminb$par

    
    if (weighting == "beta.restricted") {
      
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              w1 = 1, w2 = par["w2"],
                              theta = par["theta"],
                              m = par["m"], K = K)$tau
      
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      
    }
    if (weighting == "beta.unrestricted") {
     
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
    }
    
    
    returns <- unlist(data[y])
    
    if (gamma == TRUE) {
      # g <- c(rep(NA, times = sum(is.na((returns-par["mu"])/sqrt(tau)))),
      #        calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
      #                    alpha = par["alpha"],
      #                    beta = par["beta"],
      #                    gamma = par["gamma"],
      #                    as.numeric(na.exclude((returns-par["mu"])/sqrt(tau))),
      #                    g0 = g_zero))
      g <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         returns = as.numeric(na.exclude((returns)/sqrt(tau))),
                         g0 = g_zero))
      if(two_side==TRUE){
        e_0 <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
               calculate_shape_value(b1 = par["Pot_p1_0"], b2 = par["Pot_p2_0"],b3 = par["Pot_p3_0"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        e_1 <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
               calculate_shape_value(b1 = par["Pot_p1_1"], b2 = par["Pot_p2_1"],b3 = par["Pot_p3_1"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        e=data.frame(e_0,e_1)
        # e <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
        #        calculate_shape_value(b1 = par["Pot_p1"], b2 = par["Pot_p2"],b3 = par["Pot_p3"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        # 
        }else if(two_side==FALSE){
        e <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
               calculate_shape_value(b1 = par["Pot_p1"], b2 = par["Pot_p2"],b3 = par["Pot_p3"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
      }
    }else{
      g <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],# - par["gamma"]/2
                         alpha = par["alpha"],
                         beta = par["beta"],
                         # gamma = par["gamma"],
                         gamma = 0,
                         returns = as.numeric(na.exclude((returns)/sqrt(tau))),
                         g0 = g_zero))
      if(two_side==TRUE){
        e_0 <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
                 calculate_shape_value(b1 = par["Pot_p1_0"], b2 = par["Pot_p2_0"],b3 = par["Pot_p3_0"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        e_1 <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
                 calculate_shape_value(b1 = par["Pot_p1_1"], b2 = par["Pot_p2_1"],b3 = par["Pot_p3_1"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        e=data.frame(e_0,e_1)
      }else if(two_side==FALSE){
        e <- c(rep(NA, times = sum(is.na((returns)/sqrt(tau)))),
               calculate_shape_value(b1 = par["Pot_p1"], b2 = par["Pot_p2"],b3 = par["Pot_p3"],y = as.numeric( na.exclude(returns/sign(tau))) ) )
        }
    }
    
    if (!(var.ratio.freq %in% c("date", low.freq))) {
      if (is.null(x.two)) {
        df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau,e = e, Pot_y = Pot_y)
      }
      
    } else {
      if (is.null(x.two)) {
        df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau,e = e, Pot_y = Pot_y)
      }
      
    }
    
    # df.fitted$residuals <- unlist((df.fitted[y]-par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
    df.fitted$residuals <- unlist((df.fitted[y]) / sqrt(df.fitted$g * df.fitted$tau))
    
  }
  df.fitted$date <- as.Date(date_backup)
  
  # calculate_P_value(p.e.nlminb)
  
  # epsilon = 0.0001 * p.e.nlminb$par
  # npar=length(p.e.nlminb$par)
  # Hessian = matrix(0, ncol = npar, nrow = npar)
  # for (i in 1:npar) {
  #   for (j in 1:npar) {
  #     x1 = x2 = x3 = x4  = p.e.nlminb$par
  #     x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
  #     x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
  #     x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
  #     x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
  #     Hessian[i, j] = ( sum(lf(x1))-sum(lf(x2))-sum(lf(x3))+sum(lf(x4)) )/
  #       (4*epsilon[i]*epsilon[j])
  #   }
  # }
  # 
  # # Step 6: 创建输出结果菜单
  # se.coef = sqrt(diag(solve(Hessian)))
  # tval = p.e.nlminb$par/se.coef
  # matcoef = cbind(p.e.nlminb$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  # dimnames(matcoef) = list(names(tval), c(" Estimate",
  #                                         " Std. Error", " t value", "Pr(>|t|)"))
  # #cat("\nCoefficient(s):\n")
  # printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  # 
  
  inv_hessian <- try({
    solve(-suppressWarnings(hessian(x = par, func = function (theta) {
      -sum(lf(theta))
    }, method.args=list(eps=1e-4, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))))
  }, silent = TRUE)
  
  opg.std.err <- try({sqrt(diag(solve(crossprod(jacobian(func = function(theta) -lf(theta), x = par)))))},
                     silent = TRUE)
  if (class(opg.std.err)[1] == "try-error") {
    warning("Inverting the OPG matrix failed. No OPG standard errors calculated.")
    opg.std.err <- NA
  } else {
    opg.std.err <- opg.std.err * sqrt((mean(df.fitted$residuals^4, na.rm = TRUE) - 1) / 2)
  }
  
  if (class(inv_hessian)[1] == "try-error") {
    warning("Inverting the Hessian matrix failed. No robust standard errors calculated. Possible workaround: Multiply returns by 100.")
    rob.std.err <- NA
  } else {
    rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(func = lf, x = par)) %*% inv_hessian))
  }

  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         std.err = rob.std.err,
         broom.mgarch = data.frame(term = names(par),
                                   estimate = par,
                                   rob.std.err = rob.std.err,
                                   p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err)))),
                                   opg.std.err = opg.std.err,
                                   opg.p.value = 2 * (1 - pnorm(unlist(abs(par/opg.std.err))))),
         tau = tau,
         g = g,
         e = e,
         df.fitted = df.fitted,
         K = K,
         weighting.scheme = weighting,
         llh = p.e.nlminb$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (p.e.nlminb$value),
         # BIC = log(n) k - 2 log(L)
         y = y,
         # threshold_value = sort(ret)[floor(length(ret)*threshold_value)],
         optim = p.e.nlminb
         # , matcoef = matcoef
         )
  
  if (!is.null(x.two)) {
    output$K.two <- K.two
    output$weighting.scheme.two <- weighting.two
  }
  if (K == 0) {
    output$tau.forecast <- exp(par["m"])
  }
  
  
  # Additional output if there is a long-term component (K > 0) -------------------------------------
  if (K > 0) {
    output$variance.ratio <- 100 *
      var(log(aggregate(df.fitted$tau, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE) /
      var(log(aggregate(df.fitted$tau * df.fitted$g, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE)
    output$tau.forecast <- tau_forecast
    
    if (weighting == "beta.restricted") {
      output$est.weighting <- calculate_phi(1, w2 = par["w2"], K = K)
    }
    if (weighting == "beta.unrestricted") {
      output$est.weighting <- calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K)
    }
    if (is.null(x.two) == FALSE) {
      if (K.two > 1) {
        output$est.weighting.two <- calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two)
      }
    }
    
  }
  
  # Add class mfGARCH for employing generic functions
  class(output) <- "mfGARCH"
  output
}