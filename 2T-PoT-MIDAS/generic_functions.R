#' @export
print.mfGARCH <- function(x, ...) {
    if (class(x) != "mfGARCH") {
        stop("Obejct is not in class mfGARCH")
    } else {
        print(x$broom.mgarch)
    }
}

#' @importFrom graphics lines
#' @export
# plot.mfGARCH <- function(x, ...) {
#   if (class(x) != "mfGARCH") {
#     stop("Obejct is not in class mfGARCH")
#   }
# 
#   if (x$K == 0) {
#     plot(x = x$df.fitted["date"], y = sqrt(x$df.fitted$g),
#          type = "l",
#          xlab = colnames(x$df.fitted)[3], ylab = "vol",
#          main = "sqrt(g)", sub = "")
#   } else {
# 
#     df_plot <- aggregate(x$df.fitted, by = list(x$df.fitted[, 3]), FUN = mean)
#     #按周（也就是x$df.fitted[, 3]）分组，对每一周的所有数据去平局值
#     plot(x = df_plot[, 1], y = sqrt(df_plot$g),
#          type = "l",
#          xlab = colnames(x$df.fitted)[3], ylab = "vol",
#          main = "sqrt(tau * g) and sqrt(tau) in red", sub = "")
#     lines(x = df_plot[, 1],
#           y = sqrt(df_plot$tau),
#           col = "red")
#   }
# }

#' @export
predict.mfGARCH <- function(object, horizon = c(1:10), method = 'unnormal',fcts.tau = NULL, y.last = NULL, cond.var = NULL, cond.tau = NULL, ...) {
  if (class(object) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (is.null(cond.var)) {
    cond.var <- tail(object$g, 1)
    #tail()函数是R语言中的一个函数，用于获取向量、矩阵、表、DataFrame或函数的最后部分。
    #tail(x, n)的语法中，x是指定的数据类型，n是需要打印的行数。例如，以下代码将在名为df的数据框中获取最后3行
  }

  if (is.null(cond.tau)) {
    cond.tau <- tail(object$tau, 1)
  }

  if (is.null(fcts.tau)) {
    fcts.tau <- object$tau.forecast
  }

  if (is.null(y.last)) {
    return <- tail(object$df.fitted[object$y], 1)
  } else {
    return <- y.last
  }
  if(method == 'normal'){
    
    if (is.na(object$par["mu"])){
      ret = (return)/ sqrt(cond.tau)
      fcts.tau * as.numeric(
        1 - object$par["alpha"] - object$par["beta"] - 0/2
        + object$par["alpha"]*return^2
        + object$par["beta"]*( last(object$g) )^2)
    }
    else{if(is.na(object$par["gamma"])) {
      ret = (return - object$par["mu"])/ sqrt(cond.tau)
      fcts.tau * as.numeric(
        1 - object$par["alpha"] - object$par["beta"] - 0/2
        + object$par["alpha"]*return^2
        + object$par["beta"]*( last(object$g) )^2)
    } else {
      ret = (return - object$par["mu"])/ sqrt(cond.tau)
      fcts.tau * as.numeric(
        1 - object$par["alpha"] - object$par["beta"] - object$par["gamma"]/2
        + object$par["alpha"]*return^2
        + object$par["gamma"]*return^2*as.numeric(return<0)
        + object$par["beta"]*( last(object$g) )^2)
    }
    }
    
  }else if( method == 'unnormal' ){
  if (is.na(object$par["mu"])){
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - 0/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = 0,
                                 ret = (return)/ sqrt(cond.tau),
                                 g = cond.var))

  }
  else{if(is.na(object$par["gamma"])) {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - 0/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = 0,
                                 ret = (return - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
  } else {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - object$par["gamma"]/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = object$par["gamma"],
                                 ret = (return - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
      }
    }
  }

}

# predict.mfGARCH <- function(object, fcts.tau = NULL, y.last = NULL, cond.var = NULL, cond.tau = NULL, ...) {
#   if (class(object) != "mfGARCH") {
#     stop("Obejct is not in class mfGARCH")
#   }
#   cond.var <- tail(object$g, 1)
#     #tail()函数是R语言中的一个函数，用于获取向量、矩阵、表、DataFrame或函数的最后部分。
#     #tail(x, n)的语法中，x是指定的数据类型，n是需要打印的行数。例如，以下代码将在名为df的数据框中获取最后3行
#   cond.tau <- tail(object$tau, 1)
#   fcts.tau <- object$tau.forecast
#   return <- tail(object$df.fitted[object$y], 1)
# 
#   omega = 1 - object$par["alpha"] - object$par["beta"] - 0/2
#   alpha = object$par["alpha"]
#   beta = object$par["beta"]
#   gamma = 0
#   ret = (return)/ sqrt(cond.tau)
#   g = cond.var
# 
#   if (is.na(object$par["mu"])){
# 
#     fcts.g <- omega + alpha * ret^2 + gamma * ret^2 + beta * g
#     fcts.tau*fcts.g
#   
#   }
#   else{if(is.na(object$par["gamma"])) {
#     ret <- (return - object$par["mu"])/ sqrt(cond.tau)
#     fcts.g <- omega + alpha * ret^2 + gamma * ret^2 + beta * g
#     fcts.tau * fcts.g
# 
#   } else {
#     ret <- (return - object$par["mu"])/ sqrt(cond.tau)
#     omega = 1 - object$par["alpha"] - object$par["beta"] - object$par["gamma"]/2
#     gamma <- object$par["gamma"]
#     fcts.g <- omega + alpha * ret^2 + gamma * ret^2 + beta * g
#     fcts.tau * fcts.g
# 
#   }
#   }
#   
# }

#' This function plots the weighting scheme of an estimated GARCH-MIDAS model
#' @param x mfGARCH object obtained by fit_mfgarch
#' @importFrom graphics plot
#' @export
plot_weighting_scheme <- function(x) {
  if (class(x) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (x$weighting.scheme == "beta.restricted") {
    k = c(1:x$K)
    phi = calculate_phi(w1 = 1, w2 = x$par["w2"], K = x$K)
    plot(k, phi)

  }

  if (x$weighting.scheme == "beta.unrestricted") {
    k = c(1:x$K)
    phi = calculate_phi(w1 = x$par["w1"], w2 = x$par["w2"], K = x$K)
    plot(k, phi)
  }

}
