# JB检验
library(nortest)
library(tseries)
jb <- jarque.bera.test(modwt_df$return)
jb$p.value

# ADF检验
library(tseries)
adf <- adf.test(modwt_df$return)
adf$p.value

# install.packages("wavelets")
# install.packages("ggplot2")
# install.packages("reshape2")
summary(modwt_df$return)
sd(modwt_df$return)
library(e1071)
skewness(modwt_df$return)
library(moments)
kurtosis(modwt_df$return)


library(wavelets)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# 运行garch-midas-RV-Pot-modwt代码获取小波数据
# 子序列可视化
data=return_modwt
data$date=as.Date(modwt_df$date,format = '%Y-%m-%d')
data$return=modwt_df$return
corr(data$rv1,data$rv2,data$rv3,data$rv4,data$rv5,data$rv)
# 选择配色方案
# palette <- brewer.pal(n = 3, name = "PuBu")

# 绘制图形
p1<- ggplot(data = data, aes(x = date))+
  geom_path(aes(y = return,color='return'),size = 2,alpha = 0.7)+
  geom_path(aes(y = d1,color='d1'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = d2,color='d2'),size = 1.5,alpha = 0.7)+
  scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('return',"d1","d2")) +
  labs(x = "时间", y = "收益") +
  guides(colour = guide_legend(title = "收益序列"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1)
        ,panel.background = element_rect(fill = "white")
        ,panel.grid = element_line(color = "lightgray"))+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")
p1
p2<- ggplot(data = data, aes(x = date))+
  geom_path(aes(y = d3,color='d3'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = d4,color='d4'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = d5,color='d5'),size = 1.5,alpha = 0.7)+
  scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('d3',"d4","d5")) +
  labs(x = "时间", y = "收益") +
  guides(colour = guide_legend(title = "收益序列"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1)
        ,panel.background = element_rect(fill = "white")
        ,panel.grid = element_line(color = "lightgray"))+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")
p2

data=return_modwt
data$date=as.Date(modwt_df$date,format = '%Y-%m-%d')
data$return=modwt_df$return
# 
# # 选择配色方案
# # palette <- brewer.pal(n = 3, name = "PuBu")
# 
# # 绘制图形
# p1<- ggplot(data = data, aes(x = date))+
#   geom_path(aes(y = return,color='return'),size = 2,alpha = 0.7)+
#   geom_path(aes(y = d1,color='d1'),size = 1.5,alpha = 0.7)+
#   geom_path(aes(y = d2,color='d2'),size = 1.5,alpha = 0.7)+
#   scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('return',"d1","d2")) +
#   labs(x = "时间", y = "收益") +
#   guides(colour = guide_legend(title = "收益序列"))+
#   theme(panel.border = element_rect(fill=NA,color="black", size=1)
#         ,panel.background = element_rect(fill = "white")
#         ,panel.grid = element_line(color = "lightgray"))+
#   scale_x_date(date_breaks = "2 year", date_labels = "%Y")
# p1
# p2<- ggplot(data = data, aes(x = date))+
#   geom_path(aes(y = d3,color='d3'),size = 1.5,alpha = 0.7)+
#   geom_path(aes(y = d4,color='d4'),size = 1.5,alpha = 0.7)+
#   geom_path(aes(y = d5,color='d5'),size = 1.5,alpha = 0.7)+
#   scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('d3',"d4","d5")) +
#   labs(x = "时间", y = "收益") +
#   guides(colour = guide_legend(title = "收益序列"))+
#   theme(panel.border = element_rect(fill=NA,color="black", size=1)
#         ,panel.background = element_rect(fill = "white")
#         ,panel.grid = element_line(color = "lightgray"))+
#   scale_x_date(date_breaks = "2 year", date_labels = "%Y")
# p2

# rv可视化
data=modwt_df
data$date=as.Date(modwt_df$date,format = '%Y-%m-%d')
data$return=modwt_df$return

# 选择配色方案
# palette <- brewer.pal(n = 3, name = "PuBu")

# 绘制图形
p1<- ggplot(data = data, aes(x = date))+
  geom_path(aes(y = rv0,color='rv0'),size = 2,alpha = 0.7)+
  geom_path(aes(y = rv1,color='rv1'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = rv2,color='rv2'),size = 1.5,alpha = 0.7)+
  scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('rv0',"rv1","rv2")) +
  labs(x = "时间", y = "收益") +
  guides(colour = guide_legend(title = "收益序列"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1)
        ,panel.background = element_rect(fill = "white")
        ,panel.grid = element_line(color = "lightgray"))+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")
p1
p2<- ggplot(data = data, aes(x = date))+
  geom_path(aes(y = rv3,color='rv3'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = rv4,color='rv4'),size = 1.5,alpha = 0.7)+
  geom_path(aes(y = rv5,color='rv5'),size = 1.5,alpha = 0.7)+
  scale_colour_manual(values = c("lightgray",'#0000FF','#F53A14'),breaks=c('rv3',"rv4","rv5")) +
  labs(x = "时间", y = "收益") +
  guides(colour = guide_legend(title = "收益序列"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1)
        ,panel.background = element_rect(fill = "white")
        ,panel.grid = element_line(color = "lightgray"))+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")
p2

for (i in 0:7) {
  var_name01 <- paste0("norm_fit_model", i)  # 使用paste0函数生成变量名
  var_name02 <- paste0("t_fit_model", i)  # 使用paste0函数生成变量名
  var_name1 <- paste0("Pot_fit_model0.85_", i)  # 使用paste0函数生成变量名
  var_name2 <- paste0("Pot_fit_model0.9_", i)  # 使用paste0函数生成变量名
  var_name3<- paste0("Pot_fit_model0.95_", i)  # 使用paste0函数生成变量名
  var_name4 <- paste0("Pot_fit_model0.97_", i)  # 使用paste0函数生成变量名
  print(i)
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name01)$g*get(var_name0)$tau)[un_na])^2)))
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name02)$g*get(var_name1)$tau)[un_na])^2)))
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name1)$g*get(var_name1)$tau)[un_na])^2)))
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name2)$g*get(var_name2)$tau)[un_na])^2)))
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name3)$g*get(var_name3)$tau)[un_na])^2)))
  print(mean(abs((modwt_df$rv[un_na]-(get(var_name4)$g*get(var_name4)$tau)[un_na])^2)))
}
