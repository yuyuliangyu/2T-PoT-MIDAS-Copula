two_side_residuals_distribution <- function(x,parameters,data,threshold){
  
  # Initialization parameters
  
  sigma <- parameters[["sigma"]]
  e_0 <- parameters[["e_0"]]
  e_1 <- parameters[["e_1"]]
  k_0 <- parameters[["k_0"]]
  k_1 <- parameters[["k_1"]]
  sigma <-sqrt(sigma)
  if ( ( "threshold_0" %in% names(threshold) ) && ( "threshold_1" %in% names(threshold) ) )
  {
    threshold_0 <- threshold[["threshold_0"]]
    threshold_1 <- threshold[["threshold_1"]]
  }else{
    threshold_0 <- threshold
    threshold_1 <- threshold
  }
  
  # Calculation Threshold
  
  threshold_0 <- sort(data)[floor(length(data)*threshold_0)]
  threshold_1 <- abs(sort(data)[floor(length(data)*(1-threshold_1))])# Since the value of threshold_1 is going to be used in the code to calculate the distribution,
  # it is convenient to write it in absolute values to avoid the confusion of the -+ sign
  
  # Define the truncated empirical distribution function, that is, the empirical distribution function at less than threshold
  
  empirical_distribution <- function(x,data,threshold_0,threshold_1){
    return( sum((data<x)&(data<threshold_0)&(data>(-threshold_1)) )/sum( (data<threshold_0)&(data>(-threshold_1)) ) )
  }

  if(x>=threshold_0)
  {
    return(
      (1-k_0*(threshold_0/sigma)^(-e_0)-k_1*(threshold_1/sigma)^(-e_1))
      *empirical_distribution(x,data,threshold_0,threshold_1)
      +k_0*( (threshold_0/sigma)^(-e_0) )*( 1-(x/threshold_0)^(-e_0) )
      +k_1*( (threshold_1/sigma)^(-e_1) )
    )
  }else if( (x<threshold_0)&((-threshold_1)<x) ){
    return(
      (1-k_0*(threshold_0/sigma)^(-e_0)-k_1*(threshold_1/sigma)^(-e_1) )
      *empirical_distribution(x,data,threshold_0,threshold_1)
      +k_1*( (threshold_1/sigma)^(-e_1) )
    )
  }else if(x<=(-threshold_1)){
    return( k_1*( (threshold_1/sigma)^(-e_1) )*( 1-( 1-(-x/threshold_1)^(-e_1) ) ) )
    
  }
  
}

two_side_reverse_residuals_distribution<- function(u,parameters,data,threshold){
  
  ff <- function(theta) { (100*(two_side_residuals_distribution(theta,parameters,data,threshold)-u))^2 }
  result <-  optim(1, function(x) ff(x))$par
  return(result)
}

# Define functions to generate bilateral pot data

generate_two_side_pot_data <- function(n,parameters,data,threshold) {
  random_numbers <- n
  alpha <- parameters[["alpha"]]
  beta <- parameters[["beta"]]
  omega <- parameters[["omega"]]
  if (( "k_0" %in% names(parameters) ) ){
    k_0 <- parameters[["k_0"]]
    k_1 <- parameters[["k_1"]]
  }else{
    k_0 = 1
    k_1 = 1
  }
  
  Pot_p1_0 <- parameters[["Pot_p1_0"]]
  Pot_p2_0 <- parameters[['Pot_p2_0']]
  Pot_p3_0 <- parameters[['Pot_p3_0']]
  Pot_p1_1 <- parameters[["Pot_p1_1"]]
  Pot_p2_1 <- parameters[['Pot_p2_1']]
  Pot_p3_1 <- parameters[['Pot_p3_1']]
  threshold_0 <- threshold[["threshold_0"]]
  threshold_1 <- threshold[["threshold_1"]]
  data <- data

  e_0 <- rep( mean(calculate_shape_value(b1 = Pot_p1_0, b2 = Pot_p2_0, b3 = Pot_p3_0, (data) )), random_numbers)
  e_1 <- rep( mean(calculate_shape_value(b1 = Pot_p1_1, b2 = Pot_p2_1, b3 = Pot_p3_1, (data) )), random_numbers)

  sigma <- rep( (var(data)), random_numbers)
  residuals <- rep(NA,random_numbers) 

  random_u <- runif(1)
  
  i = 1
  probility = rep(0,random_numbers) 
  while(i < (random_numbers) ) {
    residuals[i] = two_side_reverse_residuals_distribution( random_u,c( sigma = sigma[i], e_0 = e_0[i], e_1 = e_1[i], k_0 = k_0,
                                                                        k_1 = k_1 ), data, c( threshold_0 = threshold_0, 
                                                                                              threshold_1 = threshold_1 ) )
    e_0[i+1] = exp( Pot_p1_0+Pot_p2_0*log(e_0[i])+Pot_p3_0*exp( -abs(residuals[i]) ) )
    e_1[i+1] = exp( Pot_p1_1+Pot_p2_1*log(e_1[i])+Pot_p3_1*exp( -abs(residuals[i]) ) )
    sigma[i+1] = omega+alpha*residuals[i]*residuals[i]+beta*sigma[i]
    probility[i+1]= (( ( k_0*( (sort(data)[threshold_0*length(data)]) / (sqrt(sigma[i+1]) ) )^(-e_0[i+1]) ) + k_1*( ( (sort(-1*data)[threshold_1*length(data)]) / ( sqrt(sigma[i+1]) ) )^(-e_1[i+1]) )  ) )
    
    while(  ( ( k_0*( (sort(data)[threshold_0*length(data)]) / (sqrt(sigma[i+1]) ) )^(-e_0[i+1]) ) 
              + k_1*( ( (sort(-1*data)[threshold_1*length(data)]) / ( sqrt(sigma[i+1]) ) )^(-e_1[i+1]) )  )>1){
      print("Please reset the parameters to ensure that the probability sum of the upper and lower tails is less than 1")
      print("False tail probability")
      print(( ( k_0*( (sort(data)[threshold_0*length(data)]) / (sqrt(sigma[i+1]) ) )^(-e_0[i+1]) ) + k_1*( ( (sort(-1*data)[threshold_1*length(data)]) / ( sqrt(sigma[i+1]) ) )^(-e_1[i+1]) )  ) )
      print("sigma")
      print(i)

      i<-max(i-20,0)
      print("Restart the simulation until the tail probability is less than 1 all the way through")
    }  
    random_u <- runif(1)
    i <- i + 1
  }
  residuals[random_numbers] = two_side_reverse_residuals_distribution(random_u,c( sigma = sigma[random_numbers],
                                                                                  e_0 = e_0[random_numbers],e_1 = e_1[random_numbers],
                                                                                  k_0 = k_0, k_1 = k_1 ), 
                                                                      data, c( threshold_0 = threshold_0
                                                                               , threshold_1 = threshold_1 ) )
  return(list(return = residuals, sigma = sigma, e_0 = e_0, e_1 = e_1, probility = probility ))
}

# two—sided PoT likelihood function

llh_two_side_Pot <-
  function( ret, Pot_y, k = c(k_0 = 1,k_1 = 1), omega, alpha, beta,
            threshold_value, Pot_p1, Pot_p2, Pot_p3, g_zero) {
    k_0<-k[1]
    k_1<-k[2]
    ret <- ret #-mu
    yt = sign(Pot_y)
    sigma1 <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = 0, returns = ret, g0 = 1 )
    
    sigma1 = sqrt(sigma1) 
    e_0<-calculate_shape_value(b1 = Pot_p1[1],b2 = Pot_p2[1],b3 = Pot_p3[1],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
    e_1<-calculate_shape_value(b1 = Pot_p1[2],b2 = Pot_p2[2],b3 = Pot_p3[2],(ret))#在极值理论文献中mu默认为0，但是混频考虑了mu,因此改成了(ret - mu)
    penalty0 <- ifelse( ( 1-k_0*(threshold_value[1]/sigma1)^(-e_0)-k_1*(threshold_value[2]/sigma1)^(-e_1) )>0
                        , log( 1-k_0*(threshold_value[1]/sigma1)^(-e_0)-k_1*(threshold_value[2]/sigma1)^(-e_1) ),-1e6 )
    penalty1 <- ifelse( ( ( k_0*( (threshold_value[1]) /sigma1)^(-e_0) ) + k_1*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                        , ( log(k_0)+log(e_0)+e_0*log(sigma1)-(e_0+1)*log(threshold_value[1]+Pot_y) )
                        , -1e6 )
    penalty2 <- ifelse( ( ( k_0*( (threshold_value[1]) /sigma1)^(-e_0) ) + k_1*( (threshold_value[2]) /sigma1)^(-e_1)  )<1
                        , ( log(k_1)+log(e_1)+e_1*log(sigma1)-(e_1+1)*log(threshold_value[2]+Pot_y) )
                        , -1e6 )

    -(((yt-1)/-1)*penalty0
      +yt*(ret>0)*(penalty1)
      +yt*(ret<0)*(penalty2))
  }

fit_two_side_pot_distribution <- function( ret, threshold , k_select = FALSE, sigma_select = TRUE, omega_select = TRUE ){
  Pot_y=ret
  Pot_y[ret>0] = ret[ret>0]-sort(ret)[floor(length(ret)*threshold[["threshold_0"]])]
  Pot_y[ret<0] = Pot_y[ret<0]-sort(ret)[floor(length(ret)*(1-threshold[["threshold_1"]]))]
  Pot_y[ret>0] = ifelse(Pot_y[ret>0]>0,Pot_y[ret>0],0)
  Pot_y[ret<0] = ifelse(Pot_y[ret<0]<0,abs(Pot_y[ret<0]),0)
  
  threshold_value=c( sort(ret)[floor(length(ret)*threshold[["threshold_0"]])],
                     ( -sort(ret)[floor(length(ret)*(1-threshold[["threshold_1"]]))] ) )
  if(k_select == FALSE){
    if(sigma_select == TRUE){
      if(omega_select == TRUE){
        lf <- function(p) {
          llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                           omega = p["omega"],alpha = p["alpha"], beta  = p["beta"],
                           Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                           Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                           Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]))
        }
        par.start <- c(omega = 0.01, alpha = 0.01, beta = 0.8, 
                       Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                       Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
        
        ui.opt <- rbind(c(0,-1,-1, 0, 0, 0, 0, 0, 0),
                        c(0, 1, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0,-1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 1, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 0,-1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 1))
        
        ci.opt <- c(-0.99999, 0, 0, -0.99999, 0, 0, -0.99999, 0, 0) 
      }
      else{
        lf <- function(p) {
          llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                           omega = 1-p["alpha"]-p["beta"],
                           alpha = p["alpha"], beta  = p["beta"],
                           Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                           Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                           Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]) )
        }
        par.start <- c(alpha = 0.01, beta = 0.8, 
                       Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                       Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
        
        ui.opt <- rbind(c(-1,-1, 0, 0, 0, 0, 0, 0),
                        c( 1, 0, 0, 0, 0, 0, 0, 0),
                        c( 0, 1, 0, 0, 0, 0, 0, 0),
                        c( 0, 0, 0,-1, 0, 0, 0, 0),
                        c( 0, 0, 0, 1, 0, 0, 0, 0),
                        c( 0, 0, 0, 0, 1, 0, 0, 0),
                        c( 0, 0, 0, 0, 0, 0,-1, 0),
                        c( 0, 0, 0, 0, 0, 0, 1, 0),
                        c( 0, 0, 0, 0, 0, 0, 0, 1))
        
        ci.opt <- c(-0.99999, 0, 0, -0.99999, 0, 0, -0.99999, 0, 0) 
      }
    }else{
      lf <- function(p) {
        llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                         omega = 1, alpha = 0, beta  = 0,
                         Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                         Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                         Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]))
      }
      par.start <- c(Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                     Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
      
      ui.opt <- rbind(c( 0,-1, 0, 0, 0, 0),
                      c( 0, 1, 0, 0, 0, 0),
                      c( 0, 0, 1, 0, 0, 0),
                      c( 0, 0, 0, 0,-1, 0),
                      c( 0, 0, 0, 0, 1, 0),
                      c( 0, 0, 0, 0, 0, 1))
      
      ci.opt <- c( -0.99999, 0, 0, -0.99999, 0, 0) 
    }
  }else{
    if(sigma_select == TRUE){
      if(omega_select == TRUE){
        lf <- function(p) {
          llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                           k = c(p["k_0"],p["k_1"]),
                           omega = p["omega"],alpha = p["alpha"], beta  = p["beta"],
                           Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                           Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                           Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]))
        }
        par.start <- c(k_0 = 1, k_1 = 1, omega = 0.01, alpha = 0.01, beta = 0.8, 
                       Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                       Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
        
        ui.opt <- rbind(c(0, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
        
        ci.opt <- c(-0.99999, 0, 0, -0.99999, 0, 0, -0.99999, 0, 0) 
      }
      else{
        lf <- function(p) {
          llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                           k = c(p["k_0"],p["k_1"]),
                           omega = 1-p["alpha"]-p["beta"],
                           alpha = p["alpha"], beta  = p["beta"],
                           Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                           Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                           Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]) )
        }
        par.start <- c(k_0 = 1, k_1 = 1, alpha = 0.01, beta = 0.8, 
                       Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                       Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
        
        ui.opt <- rbind(c(0, 0,-1,-1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
        
        ci.opt <- c(-0.99999, 0, 0, -0.99999, 0, 0, -0.99999, 0, 0) 
      }
    }else{
      lf <- function(p) {
        llh_two_side_Pot(ret = ret,Pot_y = Pot_y, threshold_value = threshold_value,
                         k = c(p["k_0"],p["k_1"]),
                         omega = 1, alpha = 0, beta  = 0,
                         Pot_p1 = c(p["Pot_p1_0"],p["Pot_p1_1"]),
                         Pot_p2 = c(p["Pot_p2_0"],p["Pot_p2_1"]),
                         Pot_p3 = c(p["Pot_p3_0"],p["Pot_p3_1"]))
      }
      par.start <- c(k_0 = 1, k_1 = 1, Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1,
                     Pot_p1_1 = 0.1, Pot_p2_1 = 0.8, Pot_p3_1 = 0.1)
      
      ui.opt <- rbind(c(0, 0, 0,-1, 0, 0, 0, 0),
                      c(0, 0, 0, 1, 0, 0, 0, 0),
                      c(0, 0, 0, 0, 1, 0, 0, 0),
                      c(0, 0, 0, 0, 0, 0,-1, 0),
                      c(0, 0, 0, 0, 0, 0, 1, 0),
                      c(0, 0, 0, 0, 0, 0, 0, 1))
      
      ci.opt <- c( -0.99999, 0, 0, -0.99999, 0, 0) 
    }
  }
  
  
  p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                            grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
  par <- p.e.nlminb$par
  
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
  
  p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err))))
  
  sigma1 <- calculate_g(omega = p.e.nlminb$par[['omega']], alpha = p.e.nlminb$par[['beta']]
                        , beta = p.e.nlminb$par[['alpha']], gamma = 0, returns = ret, g0 = var(ret) )
  
  
  list(par = p.e.nlminb$par, value = -p.e.nlminb$value, p.value = p.value, sigma = sigma1)
}








# One-sided

# Define the distribution function. Since pot studies tail behavior and does not know the distribution of parameters below a threshold.
# Therefore refer to zhao (2021) mentioned by Tawn (2004) semiparametric method to extend pot to residuals
one_side_residuals_distribution <- function(x,parameters,data,threshold){
  sigma <- parameters['sigma']
  e <- parameters['e']
  k <- parameters['k']
  sigma <-sqrt(sigma)
  threshold <- sort(data)[threshold*length(data)]
  
  # Define the truncated empirical distribution function, that is, the empirical distribution function at less than threshold
  empirical_distribution <- function(x,data,threshold){
    return(sum((data<x)&(data<threshold))/sum(data<threshold))
  }
  
  if(x>((1/k)^(-1/e)*threshold))
  {
    return(
      (1-k*(threshold/sigma)^(-e))*empirical_distribution(x,data,threshold)+k*( (threshold/sigma)^(-e) )*( 1-k*(x/threshold)^(-e) )
    )
  }else{
    return(
      (1-k*(threshold/sigma)^(-e))*empirical_distribution(x,data,threshold)
    )
  }
  
}


# Define one_side the inverse function
one_side_reverse_residuals_distribution<- function(u,parameters,data,threshold){
  
  ff <- function(theta) { (one_side_residuals_distribution(theta,parameters,data,threshold)-u)^2 }
  result <-  optim(1, function(x) ff(x))$par
  return(result)
}

# Define a function to generate one-sided pot data
generate_pot_data <- function(n,parameters,data,threshold) {
  random_numbers <- n
  alpha <- parameters[["alpha"]]
  beta <- parameters[["beta"]]
  omega <- parameters[["omega"]]
  pot1_0 <- parameters[["pot1_0"]]
  pot2_0 <- parameters[['pot2_0']]
  pot3_0 <- parameters[['pot3_0']]
  threshold <- threshold
  data <- data

  e <- rep(1.8,random_numbers)
  sigma <- rep(1,random_numbers)
  residuals <- rep(NA,random_numbers) 
  return <- rep(NA,random_numbers)
  random_u <- runif(1)
  for (i in 1:(random_numbers-1)) {
    residuals[i] = one_side_reverse_residuals_distribution( random_u,c( sigma[i], e[i], 1 ), XSHG_day$return,threshold )
    return[i] = residuals[i]*sqrt(sigma[i])
    e[i+1] = exp( pot1_0+pot2_0*log(e[i])+pot3_0*( -abs(residuals[i]) ) )
    sigma[i+1] = omega+alpha*sigma[i]+beta*return[i]*return[i]
    while( ( ( (sort(data)[threshold*length(data)]) /sigma[i+1])^(-e[i+1]) )>0.5 ){
      print("Please reset the parameters to ensure that the probability of up and down tails is less than 0.5")
      residuals[i] = one_side_reverse_residuals_distribution( random_u,c( sigma[i], e[i], 1 ), XSHG_day$return,threshold )
      return[i] = residuals[i]*sqrt(sigma[i])
      e[i+1] = exp( pot1_0+pot2_0*log(e[i])+pot3_0*( -abs(residuals[i]) ) )
      sigma[i+1] = omega+alpha*sigma[i]+beta*return[i]*return[i]
    }
    random_u <- runif(1)
  }
  residuals[random_numbers] = one_side_reverse_residuals_distribution(random_u,c( sigma[random_numbers], e[random_numbers], 1 ), data, threshold )
  return[random_numbers] = residuals[random_numbers]*sqrt(sigma[random_numbers])
  return(list(return, sigma, e ))
}

llh_Pot <- function(Pot_y, ret, k = 1, sigma = 1, threshold_value,
                    omega, alpha, beta, Pot_p1, Pot_p2, Pot_p3) {
  
  sigma1 = calculate_g( omega = omega,alpha = beta, beta = alpha, gamma = 0 ,returns = ret, g0 = var(ret) )
  # Here take alpha = beta, beta = alpha because it's reversed in calculate_g, and to avoid a drastic change, just swap the order to avoid errors
  sigma1 = sqrt(sigma1) 
  yt = sign(Pot_y)
  e <- calculate_shape_value(b1 = Pot_p1,b2 = Pot_p2,b3 = Pot_p3,(ret)) # In the literature on extreme value theory mu defaults to 0,
  # but mixing takes mu into account, so it is changed to (ret - mu)
  
  penalty <- ifelse( ( 1-k*(threshold_value/sigma1)^(-e) )>0, log( 1-k*(threshold_value/sigma1)^(-e) ),-1e6)
  -(((yt-1)/-1)*penalty+yt*(log(k)+log(e)+e*log(sigma1)-(e+1)*log(threshold_value+Pot_y)))#注意区分小k和大K
  
}

fit_pot_distribution <- function( ret, threshold ,omega_select){
  
  Pot_y = ret-sort(ret)[floor(length(ret)*threshold)]
  Pot_y = ifelse(Pot_y>0,Pot_y,0)
  threshold_value=sort(ret)[floor(length(ret)*threshold)]
  if(omega_select==TRUE){
    lf <- function(p) {
      llh_Pot(ret = ret,Pot_y = Pot_y, threshold_value=threshold_value,
              omega = p["omega"],alpha = p["alpha"], beta  = p["beta"], Pot_p1 = c(p["Pot_p1_0"]), Pot_p2 = c(p["Pot_p2_0"]), Pot_p3 = c(p["Pot_p3_0"]))
    }
    par.start <- c(omega = 0.01, alpha = 0.8, beta = 0.01, Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1)
    
    ui.opt <- rbind(c(0,-1,-1, 0, 0, 0),
                    c(0, 1, 0, 0, 0, 0),
                    c(0, 0, 1, 0, 0, 0),
                    c(0, 0, 0, 0,-1, 0),
                    c(0, 0, 0, 0, 1, 0),
                    c(0, 0, 0, 0,-1,-1))
    
    ci.opt <- c(-0.99999, 0.7, 0, -1, 0.6, -0.99999) 
  }else{
    lf <- function(p) {
      llh_Pot(ret = ret,Pot_y = Pot_y, threshold_value=threshold_value,
              omega = 1-p["alpha"]- p["beta"],alpha = p["alpha"],
              beta  = p["beta"], Pot_p1 = c(p["Pot_p1_0"]), Pot_p2 = c(p["Pot_p2_0"]), Pot_p3 = c(p["Pot_p3_0"]))
    }
    par.start <- c( alpha = 0.8, beta = 0.01, Pot_p1_0 = 0.1, Pot_p2_0 = 0.8, Pot_p3_0 = 0.1)
    
    ui.opt <- rbind(c(-1,-1, 0, 0, 0),
                    c( 1, 0, 0, 0, 0),
                    c( 0, 1, 0, 0, 0),
                    c( 0, 0, 0,-1, 0),
                    c( 0, 0, 0, 1, 0),
                    c( 0, 0, 0, 0, 1))
    
    ci.opt <- c(-0.99999, 0, 0, -0.99999, 0, 0) 
  }
  
  
  p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                            grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
  p.e.nlminb$par
}