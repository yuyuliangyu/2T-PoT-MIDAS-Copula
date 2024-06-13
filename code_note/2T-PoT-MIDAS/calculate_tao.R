sum_tau_fcts <- function(i, m, theta, phivar, covariate, K) {
  exponential <- m
  for (j in 1:K) {
    exponential <- exponential + theta * phivar[j] * covariate[i - j]
  }
  return(exponential)
}

sum_tau <- function(m, theta, phivar, covariate, K) {
  n <- length(covariate) - K
  exponential <- numeric(n)
  for (i in 1:n) {
    exponential[i] <- m
    for (j in 1:K) {
      exponential[i] <- exponential[i] + theta * phivar[j] * covariate[K + i - j]
      #这里为什么是K+i-j是因为K+i代表了公式里面的时间，但是考虑到了之后阶数，那么k+i不能取全部的t(1到T)，只能取k+1到T
      #因此i在这里是1到T-k
    }
  }
  return(exponential)
}