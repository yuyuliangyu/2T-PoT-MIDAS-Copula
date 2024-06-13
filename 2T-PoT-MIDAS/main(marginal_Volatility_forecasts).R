rm(list=ls())
library(foreach)
library(dplyr)
library(lubridate)
library(numDeriv)
library(Rcpp)
library(ggplot2)
library(openxlsx)
library(highfrequency)

source('helper_function.R')
source('fit_fmgarch.R')
source('generic_functions.R')
source('Pot_fit_garch.R')

################################################################################################################
# Import data

energy_car_1min<- read.csv(".\\financial_data\\energy_car.csv")
energy_car_1min<-energy_car_1min[,c(1,3)]

#rename

names(energy_car_1min)[1] <- "date"
names(energy_car_1min)[2] <- "price"

# Load the necessary libraries to process the data

library(dplyr)
library(xts)
library(lubridate)

# Converting specific trading days to date objects

date <- as.POSIXct(energy_car_1min[,'date'])
energy_car_1min[,'date']=date
energy_car_1min[2:nrow(energy_car_1min),'return'] = 100*log(energy_car_1min[2:nrow(energy_car_1min),'price'])-100*log(energy_car_1min[1:nrow(energy_car_1min)-1,'price'])

# Get 5min data

data_xts <- xts(energy_car_1min$price, order.by = energy_car_1min$date)
energy_car_5min <- to.period(data_xts, period = "minutes", k = 5)
energy_car_5min <- cbind(date = as.POSIXct(index(energy_car_5min)), as.data.frame(energy_car_5min))
energy_car_5min <- energy_car_5min[,c(1,5)]
names(energy_car_5min)[1] <- "date"
names(energy_car_5min)[2] <- "return"
energy_car_5min[2:nrow(energy_car_5min),'return'] = 100*log(energy_car_5min[2:nrow(energy_car_5min),'return'])-100*log(energy_car_5min[1:nrow(energy_car_5min)-1,'return'])
energy_car_5min = energy_car_5min[2:nrow(energy_car_5min), ]
energy_car_5min$day <- ceiling_date(energy_car_5min$date, "day") -days(1)

# Get 1day data

energy_car_day <- to.period(data_xts, period = "day", k = 1)
energy_car_day <- cbind(date = as.POSIXct(index(energy_car_day)), as.data.frame(energy_car_day))
energy_car_day <- energy_car_day[,c(1,5)]
names(energy_car_day)[1] <- "date"
names(energy_car_day)[2] <- "return"
energy_car_day[2:nrow(energy_car_day),'return'] = 100*log(energy_car_day[2:nrow(energy_car_day),'return'])-100*log(energy_car_day[1:nrow(energy_car_day)-1,'return'])
energy_car_day = energy_car_day[2:nrow(energy_car_day), ]

# 使用week()函数获取所属的周

last_day_of_week <- ceiling_date(energy_car_day$date, "week") - days(1)
last_day_of_month <- ceiling_date(energy_car_day$date, "month") - days(1)
energy_car_day[,'week']=last_day_of_week-days(1)
energy_car_day[,'month']=last_day_of_month-days(1)
df=energy_car_day
low.freq='month'

# Normalized data

df$date <- gsub("/", "-", df$date)
df$date <- format(as.Date(df$date), "%Y-%m-%d")

#############################################################################################################
#############################################################################################################
#############################################################################################################

# wavelet transform
# Extracting the low-frequency signal from the modwt decomposition, the RV at different scales will be calculated from the signal

library(waveslim)
calculate_modwt<-function(df,y='return', filter, n){
  return_modwt=modwt(t(df[y]), wf = filter, n.levels =  n, boundary = "periodic")
  return_modwt
}

modwt_df=df

rv0 <- modwt_df %>%
  group_by(modwt_df[[low.freq]]) %>%
  mutate( sum_of_squares = sum(return^2))

return_modwt = calculate_modwt(modwt_df,y='return',filter="d8", n=7)
return_modwt = as.data.frame(return_modwt[1:length(return_modwt)])
return_modwt[[low.freq]]=modwt_df[[low.freq]]

# https://blog.csdn.net/weixin_47697584/article/details/121953827
# Principles of wavelet transform reconstruction

rv1 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d1^2))
rv2 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d2^2))
rv3 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d3^2))
rv4 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d4^2))
rv5 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d5^2))
rv6 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d6^2))
rv7 <- return_modwt %>%
  group_by(return_modwt[[low.freq]]) %>%
  mutate( sum_of_squares = sum(d7^2))

modwt_df$rv0<-rv0[['sum_of_squares']]
modwt_df$rv1<-rv1[['sum_of_squares']]
modwt_df$rv2<-rv2[['sum_of_squares']]
modwt_df$rv3<-rv3[['sum_of_squares']]
modwt_df$rv4<-rv4[['sum_of_squares']]
modwt_df$rv5<-rv5[['sum_of_squares']]
modwt_df$rv6<-rv6[['sum_of_squares']]
modwt_df$rv7<-rv7[['sum_of_squares']]

# Calculating true volatility

sum_modwt_df <- energy_car_5min %>%
  group_by(energy_car_5min[['day']]) %>%
  summarise(sum_of_squares = sum(return^2))
modwt_df$rv = unlist(sum_modwt_df["sum_of_squares"])[2:length(unlist(sum_modwt_df["sum_of_squares"]))]
plot(modwt_df$rv)

# We can also replace RV with other metrics, such as RBV, MedRV, RVKernel

library(data.table)
data <- data.table(DT = energy_car_5min$date,onemin = energy_car_5min$return)
# data <- list(data)
RBV <- rBPCov(rData = data[,list(DT,onemin)],makeReturns = FALSE,makePsd = TRUE)  #RBV
MedRV<-rMedRVar(rData = data[,list(DT,onemin)], makeReturns = FALSE,makePsd = TRUE) 
RVKernel <- rKernelCov(rData = data[,list(DT,onemin)],makePsd = TRUE)
RBV <- RBV$BPV[2:nrow(RBV)]
MedRV <- MedRV$rMedRVar[2:nrow(MedRV)]
RVKernel <- RVKernel$RK[2:nrow(RVKernel)]
RBV[RBV<0] <- 0
MedRV[MedRV<0] <- 0
RVKernel[RVKernel<0] <- 0 
modwt_df$RBV <- RBV 
modwt_df$MedRV <- MedRV 
modwt_df$RVKernel <- RVKernel 
plot(RBV)
plot(MedRV)
plot(RVKernel)

# Save preprocessed data information
save(modwt_df, file = "energy_car_modwt_df.RData")
# In-sample prediction

source('helper_function.R')
source('Pot_fit_garch.R')
source('fit_fmgarch.R')

# Initialize lag order
K=24
# j=1

# Initializing the True Volatility Indicator
real_volatility <- c('rv','RBV','MedRV','RVKernel')

# Initialize the threshold set
set_threshold <- c('0.8','0.85','0.9','0.95')

# Initialize the set of loss functions
in_sample_error=data.frame()

# Initialize the model set
set_fit_model <- list()

# In-sample estimation and preservation of model data content(Here we only take i=0,
# i.e., no wavelet transform is used, and the exogenous variable is the RV calculated from the original yield series)
{
  # in-sample_prediction
  # for (i in 0:(sum(grepl( "rv", names(modwt_df)))-5) ) {
  i=0
  {
    
    var_name001 <- paste0("norm_fit_model_", i)  
    var_name002 <- paste0("t_fit_model_", i)  
    
    set_fit_model[[var_name001]] <- fit_mfgarch(data = modwt_df, y = "return", x = paste0("rv", i), low.freq = low.freq, K = K, distribution='norm')  
    set_fit_model[[var_name002]] <- fit_mfgarch(data = modwt_df, y = "return", x = paste0("rv", i), low.freq = low.freq, K = K, distribution='t')  
    # 循环拟合Pot模型
    for (s in 1:length(set_threshold)) {
      # two side
      name <- paste0('two_side_Pot_fit_model',set_threshold[s],'_',i)
      print(name)
      set_fit_model[[name]] <- Pot_fit_garch(data = modwt_df, y = "return", x = paste0("rv", i), low.freq = low.freq, K = K,
                                             threshold_value =as.numeric(set_threshold[s]),pot_k = TRUE,two_side = TRUE,gamma = TRUE)  
      print(set_fit_model[[name]][["broom.mgarch"]][["p.value"]])
      # one side
      name <- paste0('one_side_Pot_fit_model',set_threshold[s],'_',i)
      print(name)
      set_fit_model[[name]] <- Pot_fit_garch(data = modwt_df, y = "return", x = paste0("rv", i), low.freq = low.freq, K = K,
                                             threshold_value =as.numeric(set_threshold[s]),pot_k = TRUE,two_side = FALSE,gamma = TRUE)  
      print(set_fit_model[[name]][["broom.mgarch"]][["p.value"]])
    }
    
    # Cyclic Calculation of Positronic Distribution, t-Distribution Modeling Errors
    for ( w in 1:( sum( grepl( paste0("model_",i,"$"), names(set_fit_model) ) ) ) ) {
      name=names(set_fit_model)[ grepl( paste0("model_",i,"$"), names(set_fit_model) ) ][w]
      print(name)
      new_column <- calculate_error( modwt_df, set_fit_model[[name]]$g*set_fit_model[[name]]$tau)
      new_column['bic'] <- set_fit_model[[name]]$bic
      rownames(new_column) <- name 
      in_sample_error <- rbind(in_sample_error, new_column)
      # in_sample_error <- rbind(in_sample_error, as.data.frame(row.names = name,t(new_column)))
    }
    
    # Cyclic calculation of Pot error
    for ( w in 1:( sum( grepl( paste0("^one_side_Pot..*_",i,"$"), names(set_fit_model) ) ) ) ) {
      name=names(set_fit_model)[ grepl( paste0("^one_side_Pot.*_",i,"$"), names(set_fit_model) ) ][w]
      print(name)
      new_column <- calculate_error( modwt_df = modwt_df, vector = set_fit_model[[name]]$g*set_fit_model[[name]]$tau)
      new_column['bic'] <- set_fit_model[[name]]$bic
      rownames(new_column) <- name 
      in_sample_error <- rbind(in_sample_error, new_column)
    }
    for ( w in 1:( sum( grepl( paste0("^two_side_Pot.*_",i,"$"), names(set_fit_model) ) ) ) ) {
      name=names(set_fit_model)[ grepl( paste0("^two_side_Pot.*_",i,"$"), names(set_fit_model) ) ][w]
      print(name)
      new_column <- calculate_error( modwt_df = modwt_df, vector = set_fit_model[[name]]$g*set_fit_model[[name]]$tau)
      new_column['bic'] <- set_fit_model[[name]]$bic
      rownames(new_column) <- name 
      in_sample_error <- rbind(in_sample_error, new_column)
    }
  }
  
  rownames(in_sample_error[order(in_sample_error$rv_MAE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$rv_MSE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$rv_HMAE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$rv_HMSE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$rv_R2log)[1:15], ])
  
  rownames(in_sample_error[order(in_sample_error$RBV_MSE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$MedRV_MSE)[1:15], ])
  rownames(in_sample_error[order(in_sample_error$RVKernel_MSE)[1:15], ])
  
  
  write.csv(in_sample_error,'.\\in-sample_prediction\\energy_car\\in_sample_error.csv')
  save( set_fit_model, file = '.\\in-sample_prediction\\energy_car\\set_fit_model.RData')
}

# Get the parameters for all elements in the set_fit_model list, in order to save the parameters to a csv file

in_sample_estimation_result <- sapply(1:length(set_fit_model), function(i) {
  params <- set_fit_model[[i]]$par
  std_error <- set_fit_model[[i]]$broom.mgarch$rob.std.err
  p_value <- set_fit_model[[i]]$broom.mgarch$p.value
  llh <- rep(set_fit_model[[i]]$llh, length(p_value))
  bic <- rep(set_fit_model[[i]]$bic, length(p_value))
  
  data.frame(model = names(set_fit_model[i]), params = params, std_error = std_error, p_value = p_value, llh = llh, bic = bic)
})

# Save parameters to a csv file

library(openxlsx)
wb <- createWorkbook()

# Loop through each model result
for (i in 1:length(set_fit_model)) {
  model_name <- names(set_fit_model)[i]
  in_sample_estimation_result <- data.frame(
    params_name = names(set_fit_model[[i]]$par),
    params = set_fit_model[[i]]$par,
    std_error = set_fit_model[[i]]$broom.mgarch$rob.std.err,
    p_value = set_fit_model[[i]]$broom.mgarch$p.value,
    llh = rep(set_fit_model[[i]][["llh"]], length(set_fit_model[[i]]$par)),
    bic = rep(set_fit_model[[i]][["bic"]], length(set_fit_model[[i]]$par))
  )
  
  # Save the results to a different sheet in the Excel file
  addWorksheet(wb, sheetName = model_name)
  writeData(wb, sheet = model_name, x = in_sample_estimation_result)
}

# Save Excel file
saveWorkbook(wb, ".\\in-sample_prediction\\energy_car\\in_sample_estimation_result.xlsx")

############################################################################################################

# in-sample_prediction plot

{
model <- set_fit_model$norm_fit_model_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'N', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave('.\\in-sample_prediction\\energy_car\\norm_fit_model_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$t_fit_model_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'T', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\t_fit_model_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$two_side_Pot_fit_model0.8_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = '2PoT0.8', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\two_side_Pot_fit_model0.8_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$two_side_Pot_fit_model0.85_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = '2PoT0.85', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\two_side_Pot_fit_model0.85_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$two_side_Pot_fit_model0.9_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = '2PoT0.9', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\two_side_Pot_fit_model0.9_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$two_side_Pot_fit_model0.95_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = '2PoT0.95', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave('.\\in-sample_prediction\\energy_car\\two_side_Pot_fit_model0.95_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$one_side_Pot_fit_model0.8_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'PoT0.8', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave('.\\in-sample_prediction\\energy_car\\one_side_Pot_fit_model0.8_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$one_side_Pot_fit_model0.85_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'PoT0.85', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\one_side_Pot_fit_model0.85_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$one_side_Pot_fit_model0.9_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'PoT0.9', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
ggsave('.\\in-sample_prediction\\energy_car\\one_side_Pot_fit_model0.9_0(in-sample_prediction).png',plot = p)

model <- set_fit_model$one_side_Pot_fit_model0.95_0
data <- data.frame( volatility = sqrt(model$tau*model$g), tau = sqrt(model$tau), g = sqrt(model$g))
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = volatility, color = 'sqrt(tau*g)'), size = 0.7) +
  geom_line(aes(y = tau, color = 'sqrt(tau)'), size = 0.8) +
  # geom_point(aes(y = real, color = 'real'), size = 1)+
  labs(title = 'PoT0.95', x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('sqrt(tau*g)' = 'black', 'sqrt(tau)' = 'red')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave('.\\in-sample_prediction\\energy_car\\one_side_Pot_fit_model0.95_0(in-sample_prediction).png',plot = p)
}

# in-sample_prediction PoT tail index plot
{
model <- set_fit_model$one_side_Pot_fit_model0.85_0
data <- data.frame( e = model$e[!is.na(model$e)])
p <- ggplot(data, aes(x = 1:nrow(data))) +
  geom_line(aes(y = e, color = 'e'), size = 0.7) +
  labs(title = 'PoT0.85', x = "Time") +
  scale_color_manual(name = NULL,values = c('e' = 'black')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave('.\\in-sample_prediction\\energy_car\\one_side_Pot_fit_model0.85_0_dynamic_tail_index(in-sample_prediction).png',plot = p)
# 
model <- set_fit_model$two_side_Pot_fit_model0.85_0
data <- data.frame( e_0 = model$e$e_0[!is.na(model$e$e_0)], e_1 = model$e$e_1[!is.na(model$e$e_1)] )
data_e0 <- data.frame(time = 1:length(data$e_1), e_0 = data$e_0)
data_e1 <- data.frame(time = 1:length(data$e_1), e_1 = data$e_1)
plot_e0 <- ggplot(data_e0, aes(x = time, y = e_0)) +
  geom_line(color = 'black', size = 0.7) +
  labs(title = '2PoT0.85‘s e_up and e_down', x = "", y = 'e_up' ) +
  theme_minimal() + 
  theme( plot.title = element_text(hjust = 0.5) , 
         axis.title.x = element_blank(),  # Remove x-axis title
         axis.ticks.x = element_blank(),   # Remove x-axis ticks
         axis.text.x = element_blank(),   # Remove x-axis labels
         )# Remove x-axis ticks 
plot_e1 <- ggplot(data_e1, aes(x = time, y = e_1)) +
  geom_line(color = 'red', size = 0.7) +
  labs(title = '', x = "Time", y = 'e_down') +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
library(ggpubr)
p <- ggarrange(plot_e0, plot_e1, nrows = 2, ncol = 1)
ggsave('.\\in-sample_prediction\\energy_car\\two_side_Pot_fit_model0.85_0_dynamic_tail_index(in-sample_prediction).png',plot = p)
}
############################################################################################################
############################################################################################################
############################################################################################################

# Selection of models with small errors for rolling forecasts
# in_sample_error<- read.csv('C:.\\in-sample_prediction\\energy_car\\in_sample_error.csv')
# 
# order_set <- c()
# name <- names(in_sample_error)[grepl( "_", names(in_sample_error))]
# for (i in 1:sum(grepl( "_", names(in_sample_error))) ) {
#   order_set <- rbind(order_set, as.vector(in_sample_error[order(in_sample_error[[name[i]]] )[1:3],1]) )#所有误差排在前15的模型集合
# }
# order_set <- rep(order_set,1,length(order_set))
# order_set <- unique(order_set)

# Initialize the model for rolling forecasts

order_set <- as.vector( c("one_side_Pot_fit_model0.95_0","one_side_Pot_fit_model0.9_0", "one_side_Pot_fit_model0.85_0","one_side_Pot_fit_model0.8_0",
                          "two_side_Pot_fit_model0.95_0","two_side_Pot_fit_model0.9_0","two_side_Pot_fit_model0.85_0","two_side_Pot_fit_model0.8_0",
                          "norm_fit_model_0", "t_fit_model_0") )

# Initialize Rolling Prediction Parameters

df=modwt_df
T=nrow(df)
window_length=floor(T*0.9)
i=window_length
method = 'normal'

############################################################################################################

# Using multi-threaded arithmetic

library(stringr)
library(parallel)
library(progress)
#指定线程数
clus <- makeCluster(detectCores() - 10 )

clusterEvalQ(clus, {
  library(stringr)
  library(Rcpp)
  library(data.table)
}
)
clusterEvalQ(clus, source("./helper_function.R"))
clusterEvalQ(clus, source("./Pot_fit_garch.R"))
clusterEvalQ(clus, source("./generic_functions.R"))
clusterEvalQ(clus, source("./asymptotic normality_helper_function.R"))
clusterEvalQ(clus, source("./fit_fmgarch.R"))


clusterExport(clus, varlist = "df")
clusterExport(clus, varlist = "window_length")
clusterExport(clus, varlist = "low.freq")
clusterExport(clus, varlist = "K")
clusterExport(clus, varlist = "order_set")
clusterExport(clus, varlist = "T")
clusterExport(clus, varlist = "method")



result <- parLapply(clus, ( window_length ):(T-1) ,fun = rolling_predict_volatility )
stopCluster(clus)


# Extraction of forecast data

# Forecasted volatility 
predict <- as.data.frame(t((sapply( 1:length(result), function(i) { as.matrix( result[[i]][[1]] ) } )  )))
# Forecasting long-term volatility
tau_forecast <- as.data.frame(t((sapply( 1:length(result), function(i) { as.matrix( result[[i]][[2]] ) } )  )))
# Predicting the upper tail index
uptail_parameters <- as.data.frame(t((sapply( 1:length(result), function(i) { as.matrix( result[[i]][[3]] ) } )  )))
# Predicting the down tail index
downtail_parameters <- as.data.frame(t((sapply( 1:length(result), function(i) { as.matrix( result[[i]][[4]] ) } )  )))
# Predicting the tail index
tail_parameters <- as.data.frame(t((sapply( 1:length(result), function(i) { as.matrix( result[[i]][[5]] ) } )  )))

colnames(predict) <- order_set

# predict=list()
# 
# for (j in 1:length(order_set)) {
#   pre <- get(paste0('predict_',order_set[j]))
#   predict[[order_set[j]]] <- as.numeric(pre)
#   print(get(paste0('predict_',order_set[j])))
# }

# Saving data for forecasting volatility

write.csv(predict, paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_prediction(length',window_length,').csv'))
write.csv(tau_forecast, paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_tau_forecast(length',window_length,').csv'))
write.csv(uptail_parameters, paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_uptail_parameters(length',window_length,').csv'))
write.csv(downtail_parameters, paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_downtail_parameters(length',window_length,').csv'))
write.csv(tail_parameters, paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_tail_parameters(length',window_length,').csv'))

############################################################################################################
############################################################################################################
############################################################################################################

# Calculating Rolling Error

in_sample_error <- read.csv('.\\in-sample_prediction\\energy_car\\in_sample_error.csv')
in_sample_error <- in_sample_error[,2:ncol(in_sample_error)]
predict <- read.csv(paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_prediction(length',window_length,').csv'))
predict <- predict[,2:ncol(predict)]

error <- as.data.frame(matrix(nrow= ncol(predict),ncol = ncol(in_sample_error)-1 ))
rownames(error) <- names(predict)
names(error) <- names(in_sample_error)[1:(length(in_sample_error)-1)]

for (i in 1:ncol(predict)) {
  error[i,] <- calculate_error(modwt_df,as.vector(predict[,i]))
}
ncol(error)

write.csv(error,paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_prediction_error(length',window_length,').csv'))

############################################################################################################
############################################################################################################
############################################################################################################

# View a chart of predicted volatility

plot( modwt_df$rv[(window_length+1 ):T])
plot( modwt_df$RBV[(window_length+1 ):T])
plot( modwt_df$MedRV[(window_length+1 ):T])
plot( predict$norm_fit_model_0)


# Plotting real vs. predicted values
{
data <- data.frame(
  predict = sqrt(predict[["two_side_Pot_fit_model0.85_0"]]),real = sqrt(modwt_df$MedRV[(window_length+1 ):T])
)
p <-
  ggplot(data, aes(x = 1:nrow(data))) +
  geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.8) +
  geom_line(aes(y = predict, color = 'N0'), linewidth = 0.8) +
  labs(title = paste0(
    # colnames(predict)[3],
    '2PoT0.85  vs  real'), x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('real' = 'red', 'N0' = 'black')) +
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
p
ggsave(paste0('.\\out_of_sample_prediction\\energy_car\\two_side_Pot_fit_model0.85_0(length',window_length,').png'),plot = p)
# 
data <- data.frame(
  predict = sqrt(predict[["one_side_Pot_fit_model0.85_0"]]),real = sqrt(modwt_df$MedRV[(window_length+1 ):T])
)
p <-
  ggplot(data, aes(x = 1:nrow(data))) +
  geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
  geom_line(aes(y = predict, color = 'PoT0.85_5'), linewidth = 0.8) +
  labs(title = paste0(
    # colnames(predict)[3],
    'PoT0.85  vs  real'), x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('real' = 'red', 'PoT0.85_5' = 'black')) +
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
p
ggsave(paste0('.\\out_of_sample_prediction\\energy_car\\one_side_Pot_fit_model0.85_0(length',window_length,').png'),plot = p)

# 
data <- data.frame(
  predict = sqrt(predict[['norm_fit_model_0']]),real = sqrt(modwt_df$MedRV[(window_length+1 ):T])
)
p <-
  ggplot(data, aes(x = 1:nrow(data))) +
  geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
  geom_line(aes(y = predict, color = 'N6'), linewidth = 0.8) +
  labs(title = paste0(
    # colnames(predict)[3],
    'N  vs  real'), x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('real' = 'red', 'N6' = 'black')) +
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
p
ggsave(paste0('.\\out_of_sample_prediction\\energy_car\\norm_fit_model_0(length',window_length,').png'),plot = p)

# 
data <- data.frame(
  predict = sqrt(predict[["t_fit_model_0"]]),real = sqrt(modwt_df$MedRV[(window_length+1 ):T])
)
p <-
  ggplot(data, aes(x = 1:nrow(data))) +
  geom_point(aes(y = real, color = 'real'), shape = 22, size = 0.5) +
  geom_line(aes(y = predict, color = 'T0'), linewidth = 0.8) +
  labs(title = paste0(
    # colnames(predict)[3],
    'T  vs  real'), x = "Time", y = "Volatility") +
  scale_color_manual(name = NULL,values = c('real' = 'red', 'T0' = 'black')) +
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
p
ggsave(paste0('.\\out_of_sample_prediction\\energy_car\\t_fit_model_0(length',window_length,').png'),plot = p)
}


# MCS predictive capability test
library(MCS)
library(readxl)

predict <- read.csv(paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_prediction(length',window_length,').csv'))
predict <- as.data.frame(predict)
predict <- predict[,2:ncol(predict)]

perform_MCS <- function(predict, real_volatility, alpha = 0.2, B = 5000, statistic = 'Tmax', cl = NULL) {
  
  error_matrix_MAE <- as.data.frame(sapply(1:ncol(predict), function(x) abs(as.vector(predict[,x]) - real_volatility )))
  error_matrix_MSE <- error_matrix_MAE^2
  error_matrix_HMAE <- as.data.frame(sapply(1:ncol(predict),function(x) abs( 1 - as.vector( (predict[,x]) ) / real_volatility ) ) )
  error_matrix_HMSE <- error_matrix_HMAE^2
  error_matrix_Rlog <- as.data.frame(sapply(1:ncol(predict),function(x) ( log( as.vector( (predict[,x]) ) / real_volatility ) )^2 ) )
  colnames(error_matrix_MAE) <- names(predict)
  colnames(error_matrix_MSE) <- names(predict)
  colnames(error_matrix_HMAE) <- names(predict)
  colnames(error_matrix_HMSE) <- names(predict)
  colnames(error_matrix_Rlog) <- names(predict)
  
  MCSprocedure_helper <- function(Loss, alpha, B, statistic, cl) {
    MCSprocedure(Loss = Loss, alpha = alpha, B = B, statistic = statistic, cl = cl)
  }
  print('MAE')
  MCS_MAE <- MCSprocedure_helper(Loss = error_matrix_MAE, alpha = alpha, B = B, statistic = statistic, cl = cl)
  # write_xlsx(list(MCS_MAE = MCS_MAE), path = "C:\\Users\\yuyu\\Desktop\\garch-midas\\in-sample_prediction\\energy_car\MCS_results.xlsx", col_names = TRUE, format_headers = TRUE)
  print('MSE')
  MCS_MSE <- MCSprocedure_helper(Loss = error_matrix_MSE, alpha = alpha, B = B, statistic = statistic, cl = cl)
  print('HMAE')
  MCS_HMAE <- MCSprocedure_helper(Loss = error_matrix_HMAE, alpha = alpha, B = B, statistic = statistic, cl = cl)
  print('HMSE')
  MCS_HMSE <- MCSprocedure_helper(Loss = error_matrix_HMSE, alpha = alpha, B = B, statistic = statistic, cl = cl)
  print('Rlog')
  MCS_Rlog <- MCSprocedure_helper(Loss = error_matrix_Rlog, alpha = alpha, B = B, statistic = statistic, cl = cl)
  
  # Save results to CSV
  return(list(MCS_MAE = capture.output(MCS_MAE), MCS_MSE = capture.output(MCS_MSE), MCS_HMAE = capture.output(MCS_HMAE), MCS_HMSE = capture.output(MCS_HMSE)
              , MCS_Rlog = capture.output(MCS_Rlog) ) ) 
}
MCS_result <- perform_MCS(predict, modwt_df$rv[(window_length + 1):T])
MCS_result <- as.data.frame(MCS_result)

predict <- read.csv(paste0('.\\out_of_sample_prediction\\energy_car\\energy_car_rolling_prediction(length',window_length,').csv'))
write.csv(MCS_result, file = paste0('.\\out_of_sample_prediction\\energy_car\\MCS_results(length',window_length,').csv'))

wb <- createWorkbook()

# Adding data with different column names to different sheet tables

for ( i in 1:length(MCS_result)) {
  addWorksheet(wb, names(MCS_result)[i] )
  writeData(wb, names(MCS_result)[i], MCS_result[[i]])
}

# save MSC test result
saveWorkbook(wb, ".\\out_of_sample_prediction\\energy_car\\MCS_result.xlsx")

# A garch model with the smallest error is selected as the final marginal model

order_set <- c()
name <- colnames(error)
error <- cbind(rownames(error),error)
for (i in 1:sum(grepl( "_", names(error))) ) {
  order_set <- rbind(order_set, as.vector(error[order(error[[name[i]]] )[1], 1]) )# Set of all models with top 15 errors
}
order_set <- rep(order_set,1,length(order_set))

# Counting of top rankings

best_model <- names(table(order_set))[which.max(table(order_set))]
print(max(table(order_set)))
print(best_model)

library(stringr)

numbers <- str_extract_all(best_model, "\\d+\\.?\\d*")[[1]]
if(grepl('Pot',best_model)){
  best_fit_model = Pot_fit_garch(data = modwt_df, y = "return", x = paste0("rv", as.numeric(numbers[2])),
                                 low.freq = low.freq, K = K,threshold_value = as.numeric(numbers[1])
                                 , pot_k = TRUE, two_side = TRUE, gamma = TRUE ) 
  
}else{
  best_fit_model=fit_mfgarch(data = modwt_df, y = "return", x = paste0("rv", as.numeric(numbers[1]))
                             , low.freq = low.freq, K = K,distribution = str_extract(best_model, "^[a-z]+" ))
}


best_fit_model$bic
best_fit_model[["broom.mgarch"]][["p.value"]]
plot(best_fit_model)
energy_car=list(data = modwt_df,model = best_model )
current_dir <- getwd()
parent_folder <- dirname(current_dir)
filename <- paste0(parent_folder,'/2T-PoT-GAS-Copula/energy_car.RData')
save(energy_car, file = filename)
load("energy_car.RData")
