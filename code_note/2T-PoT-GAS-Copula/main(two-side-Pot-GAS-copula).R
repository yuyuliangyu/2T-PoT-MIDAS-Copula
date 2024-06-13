rm(list=ls())

# Get file location
current_dir <- getwd()
parent_folder <- dirname(current_dir)
path <- paste0(parent_folder,'/2T-PoT-MIDAS')
setwd(path)

library(numDeriv)
library(Rcpp)
library(stringr)
library(foreach)
library(dplyr)
library(lubridate)
library(highfrequency)
library(xts)
library(copula)
library(VC2copula)
library(VineCopula)

# Import marginally relevant function files 
source('helper_function.R')
source('fit_fmgarch.R')
source('generic_functions.R')
source('Pot_fit_garch.R')
source('asymptotic normality_helper_function.R')

setwd(current_dir)
source('copula_helper_funciton.R')

# Initialize model parameters
K=24
low.freq='month'

method = 1 # Improved methodology
# method = 2 # original methodology

# Selection of the distribution and the corresponding data_tail 

# distribution ='2PoT'
# distribution ='PoT'
# data_tail = 'two_side'
# data_tail = 'down'
# data_tail = 'up'
model_set = list( distribution = c('2PoT', 'PoT', 'PoT') , data_tail = c('two_side', 'down', 'up') )

for( m in 1:length(model_set$distribution)){
  distribution = model_set$distribution[m]
  data_tail = model_set$data_tail[m]

  # Creating a variable storage list
  variable_list=list()
  load('energy_car.RData')
  energy_car$model <- "Pot_fit_model0.85_0"
  variable_list[['energy_car']] <- energy_car
  load("energy.RData")
  energy$model <- "Pot_fit_model0.85_0"
  variable_list[['energy']] <- energy
  
if( (distribution == 'PoT')&&( data_tail == 'down' ) ){
  variable_list[["energy_car"]][["data"]][["return"]] <- -(energy_car$data$return)
  variable_list[["energy"]][["data"]][["return"]] <- -(energy$data$return)
}else if( (distribution == 'PoT')&&( data_tail == 'up' ) ){
  variable_list[["energy_car"]][["data"]][["return"]] <- (energy_car$data$return)
  variable_list[["energy"]][["data"]][["return"]] <- (energy$data$return)
}else if( (distribution == '2PoT') ){
  variable_list[["energy_car"]][["data"]][["return"]] <- (energy_car$data$return)
  variable_list[["energy"]][["data"]][["return"]] <- (energy$data$return)
}

# Data preprocessing
date = c( variable_list[[1]][[1]]$date[1], variable_list[[1]][[1]]$date[length(variable_list[[1]][[1]]$date)] )

# Make sure both are the same length
for (i in 2:length(variable_list) ) {
  list <- variable_list[[i]][[1]]
  variable_list[[i]][[1]] <- list[ which( list == date[1] ): which(list == date[2]), ]
}

# In-sample fitting to extract suprathresholds
{
  if( distribution == '2PoT'){
    for (i in 1:length(variable_list)){
      
      numbers <- str_extract_all(variable_list[[i]][[2]], "\\d+\\.?\\d*")[[1]]
      assign( paste0( 'best_fit_garch_', i), Pot_fit_garch(data = variable_list[[i]][[1]], y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                                           low.freq = low.freq, K = K, threshold_value = as.numeric(numbers[1])
                                                           , pot_k = TRUE, two_side = TRUE, gamma = TRUE ) )
      
      ret <- get( paste0( 'best_fit_garch_', i) )$df.fitted$return
      # Calculate the marginal Pot_y
      assign(
        paste0( 'marginal_Pot_y_', i), calculate_marginal_label( get( paste0( 'best_fit_garch_', i) ),
                                                                 threshold = c( threshold_0 = as.numeric(numbers[1]), threshold_1 = as.numeric(numbers[1]) ) )[["Pot_y"]]
      )
      # Calculate the region corresponding to Pot_y,get the < negative threshold,> threshold and the label of the middle part.
      assign(
        paste0( 'marginal_Pot_y_label_', i), calculate_marginal_label( get( paste0( 'best_fit_garch_', i) ),
                                                                       threshold = c( threshold_0 = as.numeric(numbers[1]), threshold_1 = as.numeric(numbers[1]) ) )[["label"]]
      )
      # Compute the cumulative distribution variable for the marginal Pot_y
      sigma <- get( paste0( 'best_fit_garch_', i) )$df.fitted$tau * get( paste0( 'best_fit_garch_', i) )$df.fitted$g
      e_0 <- get( paste0( 'best_fit_garch_', i) )$e$e_0
      e_1 <- get( paste0( 'best_fit_garch_', i) )$e$e_1
      if( 'k_0' %in% names(get( paste0( 'best_fit_garch_', i) )$par)){
        k_0 <- get( paste0( 'best_fit_garch_', i) )$par['k_0']
        k_1 <- get( paste0( 'best_fit_garch_', i) )$par['k_1']
        parameters <- list( k_0 = k_0, k_1 = k_1, e_0 = e_0, e_1 = e_1, sigma = sigma)
      }else{
        k_0 <- 1
        k_1 <- 1
        parameters <- list( k_0 = k_0, k_1 = k_1, e_0 = e_0, e_1 = e_1, sigma = sigma)
      }
      assign(
        paste0( 'marginal_Pot_y_cdf_', i), Pot_y_cdf( ret,parameters, threshold = c( threshold_0 = as.numeric(numbers[1]),
                                                                                     threshold_1 = as.numeric(numbers[1]) ) )
      )
      
      # Calculate the residual et
      residuals <- get( paste0( 'best_fit_garch_', i) )$df.fitted$residuals
      residuals <- residuals[ !is.na(residuals) ]
      assign( paste0("marginal_residuals_", i), residuals)
      
      # Calculate the cumulative distribution corresponding to the rate of return rt
      assign( paste0("marginal_residuals_cdf_", i), unlist( sapply( 1:length( residuals ),
                                                                    
                                                                    function(k){
                                                                      NA_indice <- length(ret) - length(residuals)
                                                                      return(
                                                                        residuals_distribution( ret[ k + NA_indice ],
                                                                                                c( sigma = sigma[k+NA_indice], e_0 = e_0[k+NA_indice], e_1 = e_1[k+NA_indice]
                                                                                                   ,k_0 = as.numeric(k_0)
                                                                                                   , k_1 = as.numeric(k_1))
                                                                                                , ret, c( threshold_0 = as.numeric(numbers[1]), threshold_1 = as.numeric(numbers[1])))
                                                                      )
                                                                    } ) ) )
    }
  }else if ( distribution == 'PoT'){
    for (i in 1:length(variable_list)){
      
      
      numbers <- str_extract_all(variable_list[[i]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      assign( paste0( 'best_fit_garch_', i), Pot_fit_garch(data = variable_list[[i]][[1]], y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                                      low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                                      , pot_k = TRUE, two_side = FALSE, gamma = TRUE ))
      ret <- get( paste0( 'best_fit_garch_', i) )$df.fitted$return
      # Calculate the marginal Pot_y
      assign(
        paste0( 'marginal_Pot_y_', i), calculate_marginal_label( get( paste0( 'best_fit_garch_', i) ),
                                                                     threshold = as.numeric(numbers[1]) )[["Pot_y"]]  )   
      
      assign(
        paste0( 'marginal_Pot_y_label_', i), calculate_marginal_label( get( paste0( 'best_fit_garch_', i) ),
                                                                           threshold = c( threshold = as.numeric(numbers[1]) ) )[["label"]]
      )
      # Compute the cumulative distribution variable for the marginal Pot_y
      sigma <- get( paste0( 'best_fit_garch_', i) )$df.fitted$tau * get( paste0( 'best_fit_garch_', i) )$df.fitted$g
      e <- get( paste0( 'best_fit_garch_', i) )$e
      # parameters <- list( k = 1, e = e, sigma = sigma)
      if( 'k' %in% names(get( paste0( 'best_fit_garch_', i) )$par)){
        k <- get( paste0( 'best_fit_garch_', i) )$par['k']
        parameters <- list( k = k, e = e, sigma = sigma)
      }else{
        k <- 1
        parameters <- list( k = k, e = e, sigma = sigma)
      }
      assign(
        paste0( 'marginal_Pot_y_cdf_', i), Pot_y_cdf( ret,parameters, threshold = c( threshold = as.numeric(numbers[1]) ) )
      )

    }
  }else if ( distribution == 'norm'){
    for (i in 1:length(variable_list)){
  
      numbers <- str_extract_all(variable_list[[i]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      assign( paste0( 'best_fit_garch_', i), fit_mfgarch(data = variable_list[[i]][[1]], y = "return", x = paste0("rv", as.numeric(numbers[2])), low.freq = low.freq, K = K, distribution='norm') )
      assign( paste0( 'marginal_y_cdf_', i), na.omit(pnorm(get( paste0( 'best_fit_garch_', i) )$df.fitted$residuals ) ) )
      
    }
  }else if ( distribution == 't'){
    for (i in 1:length(variable_list)){
      numbers <- str_extract_all(variable_list[[i]][[2]], "\\d+\\.?\\d*")[[1]]# 这个表达式的主要目的是从字符串 str 中提取所有的数字，包括整数和浮点数。
      assign( paste0( 'best_fit_garch_', i), fit_mfgarch(data = variable_list[[i]][[1]], y = "return", x = paste0("rv", as.numeric(numbers[2])), low.freq = low.freq, K = K, distribution='t') )
      assign( paste0( 'marginal_y_cdf_', i), na.omit(pt( get( paste0( 'best_fit_garch_', i) )$df.fitted$residuals, df = get( paste0( 'best_fit_garch_', i) )$par[['V']] ) ) )
      
      }
  }

{
  plot(marginal_residuals_1)
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_marginal_residuals_1.png"))
  dev.off()
  plot(marginal_residuals_2)
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_marginal_residuals_2.png"))
  dev.off()
  plot(marginal_residuals_cdf_1)
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_marginal_residuals_cdf_1.png"))
  dev.off()
  plot(marginal_residuals_cdf_2)
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_marginal_residuals_cdf_2.png"))
  dev.off()
  plot(marginal_Pot_y_1, ylab = '', xlab = "time", main = "yt1 series")
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_yt1 series.png"))
  dev.off()
  plot(marginal_Pot_y_2, ylab = '', xlab = "time", main = "yt2 series")
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_yt2 series.png"))
  dev.off()
  plot(marginal_Pot_y_cdf_1, ylab = '', xlab = "time", main = "F(yt1) series")
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_F(yt1) series.png"))
  dev.off()
  plot(marginal_Pot_y_cdf_2, ylab = '', xlab = "time", main = "F(yt2) series")
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_F(yt2) series.png"))
  dev.off()
  # 绘制边际变量的散点图
  plot(marginal_Pot_y_cdf_1,marginal_Pot_y_cdf_2, ylab = '', xlab = " ", main = "F(yt1) and F(yt2) cdf scatter plot")
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_F(yt1) and F(yt2) cdf scatter plot.png"))
  dev.off()
}
  
  
  marginal_Pot_y_label <- paste0( marginal_Pot_y_label_1, marginal_Pot_y_label_2)
  marginal_residuals_label <- rep( '11', length(marginal_residuals_1))
  aa=fit_tvc(list( marginal_Pot_y_cdf_1, marginal_Pot_y_cdf_2), marginal_Pot_y_label)
  
  aa$par
  aa$broom.mgarch$p.value
  aa$broom.mgarch$rob.std.err
  aa$llh
  plot(aa$at, type = "l", ylab = '', xlab = "time", main = paste(distribution,"time-varying parameters at"))
  dev.copy(png, paste0(".\\picture\\",distribution,"_",data_tail,"_time-varying parameters at.png"))
  dev.off()
}

# Scroll calculation var
# Initialize roll parameters
df=variable_list[[1]][[1]]
T=nrow(df)
# window_length=floor(T*0.90)
window_length=floor(T*0.80)
# Initialize VaR, calculate upper tail VaR and lower tail VaR separately
VaR_level <- c( 0.9, 0.95, 0.975, 0.99, 0.995)
VaR_set <- as.vector( matrix(NA,length(VaR_level) ) )
VaR_set <- list( total_VaR_set_up = VaR_set, total_VaR_set_down = VaR_set,
                 X1_VaR_set_up = VaR_set, X1_VaR_set_down = VaR_set, X1_CoVaR_set_up = VaR_set, X1_CoVaR_set_down = VaR_set, X1_g_CoVaR_set_up = VaR_set, X1_g_CoVaR_set_down = VaR_set,
                 X2_VaR_set_up = VaR_set, X2_VaR_set_down = VaR_set, X2_CoVaR_set_up = VaR_set, X2_CoVaR_set_down = VaR_set, X2_g_CoVaR_set_up = VaR_set, X2_g_CoVaR_set_down = VaR_set)

{
# i=window_length


library(parallel)
library(progress)
  
clus <- makeCluster(detectCores() - 10 )

clusterEvalQ(clus, {
  library(stringr)
  library(Rcpp)
  library(copula)
  library(VC2copula)
  library(VineCopula)
}
)

clusterExport(clus, varlist = "variable_list")
clusterExport(clus, varlist = "window_length")
clusterExport(clus, varlist = "distribution")
clusterExport(clus, varlist = "method")
clusterExport(clus, varlist = "low.freq")
clusterExport(clus, varlist = "K")
clusterExport(clus, varlist = "VaR_level")
clusterExport(clus, varlist = "VaR_set")
clusterExport(clus, varlist = "path")
clusterExport(clus, varlist = "current_dir")

clusterEvalQ(clus, source(paste0(path,'/helper_function.R')))
clusterEvalQ(clus, source(paste0(path,'/Pot_fit_garch.R')))
clusterEvalQ(clus, source(paste0(path,'/generic_functions.R')))
clusterEvalQ(clus, source(paste0(path,'/asymptotic normality_helper_function.R')))
clusterEvalQ(clus, source(paste0(current_dir,'/copula_helper_funciton.R')))

result <- parLapply(clus, ( window_length ):(T-1) ,fun = rolling_prediction )

stopCluster(clus)

save(result, file = paste0(current_dir, "//Rolling_forecast_var//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                            ,ifelse(distribution=='PoT',data_tail,''),"(length",toString(window_length),").RData") )
}

library(segMGarch)
library(GAS)

# Get VaR of the run result
{
total_VaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$total_VaR_set_up ) } )  ) )
total_VaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$total_VaR_set_down ) } ) ) )
X1_VaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_VaR_set_up ) } )  ) )
X2_VaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_VaR_set_up ) } )  ) )
X1_VaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_VaR_set_down ) } )  ) )
X2_VaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_VaR_set_down ) } )  ) )

X1_CoVaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_CoVaR_set_up ) } )  ) )
X1_CoVaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_CoVaR_set_down ) } )  ) )
X2_CoVaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_CoVaR_set_up ) } )  ) )
X2_CoVaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_CoVaR_set_down ) } )  ) )

X1_g_CoVaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_g_CoVaR_set_up ) } )  ) )
X1_g_CoVaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X1_g_CoVaR_set_down ) } )  ) )
X2_g_CoVaR_set_up <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_g_CoVaR_set_up ) } )  ) )
X2_g_CoVaR_set_down <- as.matrix(t(sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]]$X2_g_CoVaR_set_down ) } )  ) )

at <- t(sapply( 1:length(result), function(i) { as.numeric( result[[i]][[2]] ) }))

total_return <- 0.5*variable_list[[1]][[1]]$return+0.5*variable_list[[2]][[1]]$return
total_return <- total_return[( window_length + 1 ):T]
return_1 <- variable_list[[1]][[1]]$return[( window_length + 1 ):T]
return_2 <- variable_list[[2]][[1]]$return[( window_length + 1 ):T]
}

# Plotting the forecast var

library(ggplot2)
{
total_VaR <- BacktestVaR_result(total_return, total_VaR_set_up, total_VaR_set_down)
# p <- create_plot(total_VaR_set_up, total_VaR_set_down, total_return, method, distribution, data_tail, VaR = '总体VaR')
p <- create_plot(total_VaR_set_up, total_VaR_set_down, total_return, method, distribution, data_tail, VaR = 'VaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_',
                ifelse(distribution=='PoT',data_tail,''),"forcast_tatal_VaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 ) 


# x1 VaR
x1_VaR <-BacktestVaR_result(return_1, X1_VaR_set_up, X1_VaR_set_down)
p <- create_plot(X1_VaR_set_up, X1_VaR_set_down, return_1, method, distribution, data_tail, VaR = 'x1_VaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x1_VaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 ) 

# x2 VaR
x2_VaR <-BacktestVaR_result(return_2, X2_VaR_set_up, X2_VaR_set_down)
p <- create_plot(X2_VaR_set_up, X2_VaR_set_down, return_2, method, distribution, data_tail, VaR = 'x2_VaR')
ggsave( paste0( current_dir,"//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x2_VaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 ) 

# backtesting
# x1
x1_CoVaR <-BacktestCoVaR_result(return_1, return_2, X1_CoVaR_set_up, X2_VaR_set_up, X1_CoVaR_set_down, X2_VaR_set_down)
p <- create_plot(X1_CoVaR_set_up, X1_CoVaR_set_down, return_1, method, distribution, data_tail, VaR = 'x1_CoVaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x1_CoVaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 )

# x2
x2_CoVaR <-BacktestCoVaR_result(return_2, return_1, X2_CoVaR_set_up, X1_VaR_set_up, X2_CoVaR_set_down, X1_VaR_set_down)
p <- create_plot(X2_CoVaR_set_up, X2_CoVaR_set_down, return_2, method, distribution, data_tail, VaR = 'x2_CoVaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x2_CoVaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 )
# x1
x1_g_CoVaR <-BacktestCoVaR_result(return_1, return_2, X1_g_CoVaR_set_up, X2_VaR_set_up, X1_g_CoVaR_set_down, X2_VaR_set_down)
p <- create_plot(X1_g_CoVaR_set_up, X1_g_CoVaR_set_down, return_1, method, distribution, data_tail, VaR = 'x1_g_CoVaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x1_g_CoVaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 ) 

# x2
x2_g_CoVaR <-BacktestCoVaR_result(return_2, return_1, X2_g_CoVaR_set_up, X1_VaR_set_up, X2_g_CoVaR_set_down, X1_VaR_set_down)
p <- create_plot(X2_g_CoVaR_set_up, X2_g_CoVaR_set_down, return_2, method, distribution, data_tail, VaR = 'x2_g_CoVaR')
ggsave( paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                ,ifelse(distribution=='PoT',data_tail,''),"forcast_x2_g_CoVaR","(length",toString(window_length),").png"), plot = p, width = 10, height = 8, dpi = 300 ) 
}
merged_results <- cbind( total_VaR = total_VaR, x1_VaR = x1_VaR, x2_VaR = x2_VaR, x1_CoVaR = x1_CoVaR, x2_CoVaR = x2_CoVaR, x1_g_CoVaR = x1_g_CoVaR, x2_g_CoVaR = x2_g_CoVaR )
write.csv(merged_results, paste0( current_dir, "//picture//",ifelse(method=='2','tradition','revised'),"_method","_",distribution,'_'
                                  ,ifelse(distribution=='PoT',data_tail,''),"(length",toString(window_length),").csv") )
}
