#------------------------------------#
#            Homework 2              #
#          PM2.5 Forecast            #
#                                    #
#         Mehak Shamdasani           #
#------------------------------------#
install.packages(
  c("ggfortify", "changepoint",
    "strucchange", "ggpmisc")
)
library(haven)
library(forecast)
library(fma)
library(tseries)
library(expsmooth)
library(lmtest)
library(zoo)
library(ggplot2)
library(dplyr)
library(lubridate)


# Saving File Locations and Uploading CSV File #
air_quality <- read.csv("/Users/mehak/Desktop/MSA/FALL2020/TimeSeries/Homework/PM_2_5_Raleigh2.csv")

#Prepare and Aggregate Data---------------------------------------------------------------------------
#format dates
air_quality <- air_quality %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

#aggregate PM2.5 concetration to monthly
mon_aq<- air_quality %>%
  group_by(yr = year(Date),mon = month(Date)) %>% 
  summarise(Monthly.Mean.PM2.5.Concentration = mean(Daily.Mean.PM2.5.Concentration))

#Create TS Object---------------------------------------------------------------------------
monthly_pm2 <- ts(mon_aq$Monthly.Mean.PM2.5.Concentration, start = 2014, frequency =12)

#Decompose TS Object---------------------------------------------------------------------------
decomp_stl <- stl(monthly_pm2, s.window = 7)
plot(decomp_stl)

#If anoyone wants to see the classical decomposition
a<- decompose(monthly_pm2, type = c("additive"), filter = NULL)
m<-  decompose(monthly_pm2, type = c("multiplicative"), filter = NULL)
plot(m)
plot(a)

#Create Data Split---------------------------------------------------------------------------
pm2_train= subset(monthly_pm2, end=length(monthly_pm2)-6)
pm2_valid=subset(monthly_pm2,start=length(monthly_pm2)-5)

#Decompose Training TS Object---------------------------------------------------------------------------
decomp_stl <- stl(pm2_train, s.window = 7)
plot(decomp_stl)

#Plot Trend vs Actual TS Data---------------------------------------------------------------------------
trend_pm2 <- decomp_stl$time.series[,2]
plot(pm2_train)
lines(trend_pm2, col="red")

#Plot Seasonally Adjusted vs Actual TS Data---------------------------------------------------------------------------
seas_pm2<- pm2_train - decomp_stl$time.series[,1]
plot(pm2_train)
lines(seas_pm2, col="red")

# Building a Single Exponential Smoothing Model---------------------------------------------------------------------------
SES.pm2 <- ses(pm2_train, initial = "optimal", h = 24)
summary(SES.pm2)

plot(SES.pm2, main = "Monthly PM2.5 Concentration with Simple ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
abline(v = 1992, col = "red", lty = "dashed")
round(accuracy(SES.pm2),2)

autoplot(SES.pm2)+
  autolayer(fitted(SES.pm2),series="Fitted")+ylab("Monthly PM2.5 Concentration with Simple ESM Forecast")

#Calculate MAPE on Validation
test.ses=forecast(SES.pm2,h=6)
error=pm2_valid-test.ses$mean
MAE=mean(abs(error))
MAPE=mean(abs(error)/abs(pm2_valid))
#MAE: 1.888392
#MAPE: 0.1810493

#Calculate RMSE on Validation
error_sq=(pm2_valid-test.ses$mean)^2
MSE = mean(error_sq)
RMSE = sqrt(MSE)
#RMSE: 2.287174


# Building a Linear Exponential Smoothing Model---------------------------------------------------------------------------
LES.pm2 <- holt(pm2_train, initial = "optimal", h = 24)
summary(LES.pm2)

plot(LES.pm2, main = "Monthly PM2.5 Concentration with Linear ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 1992, col = "red", lty = "dashed")

autoplot(LES.pm2)+
  autolayer(fitted(LES.pm2),series="Fitted")+ylab("Monthly PM2.5 Concentration with Holt ESM Forecast")

#Calculate MAPE on Validation
test.les=forecast(LES.pm2,h=6)
error=pm2_valid-test.les$mean
MAE=mean(abs(error))
MAPE=mean(abs(error)/abs(pm2_valid))
#MAE: 2.200177
#MAPE: 0.2074771

#Calculate RMSE on Validation
error_sq=(pm2_valid-test.les$mean)^2
MSE = mean(error_sq)
RMSE = sqrt(MSE)
#RMSE:  2.601203

#Build a Linear Damped---------------------------------------------------------------------------
LDES.pm2 <- holt(pm2_train, initial = "optimal", h = 24, damped = TRUE)
summary(LDES.pm2)

plot(LDES.pm2, main = "Monthly PM2.5 Concentration with Linear Damped ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 1992, col = "red", lty = "dashed")

autoplot(LDES.pm2)+
  autolayer(fitted(LDES.pm2),series="Fitted")+ylab("PM2.5 Concentration")

#Calculate MAPE on Validation
test.les_d=forecast(LDES.pm2,h=6)
error=pm2_valid-test.les_d$mean
MAE=mean(abs(error))
MAPE=mean(abs(error)/abs(pm2_valid))
#MAE: 2.008073
#MAPE: 0.1906637

#Calculate RMSE on Validation
error_sq=(pm2_valid-test.les_d$mean)^2
MSE = mean(error_sq)
RMSE = sqrt(MSE)
#2.405117


#Holt-Winters Exponential Smoothing Model---------------------------------------------------------------------------
#Additive---------------------------------------------------------------------------
HWES.pm2_a <- hw(pm2_train, seasonal = "additive", initial='optimal')
summary(HWES.pm2_a)

plot(HWES.pm2_a, main = "Monthly PM2.5 Concentration with Holt-Winters ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 2008.25, col = "red", lty = "dashed")

autoplot(HWES.pm2_a)+
  autolayer(fitted(HWES.pm2_a),series="Fitted")+ylab("PM2.5 Concentration")

#Calculate MAPE on Validation
test.hw_a=forecast(HWES.pm2_a,h=6)
error=pm2_valid-test.hw_a$mean
MAE=mean(abs(error))
MAPE=mean(abs(error)/abs(pm2_valid))
#MAE: 2.2477
#MAPE = 0.2285497

#Calculate RMSE on Validation
error_sq=(pm2_valid-test.hw_a$mean)^2
MSE = mean(error_sq)
RMSE = sqrt(MSE)
#RMSE: 2.324086

#Multiplicative---------------------------------------------------------------------------
HWES.pm2_m <- hw(pm2_train, seasonal = "multiplicative", initial='optimal')
summary(HWES.pm2_m)

plot(HWES.pm2_m, main = "Monthly PM2.5 Concentration with Holt-Winters ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 2008.25, col = "red", lty = "dashed")

autoplot(HWES.pm2_m)+
  autolayer(fitted(HWES.pm2_m),series="Fitted")+ylab("PM2.5 Concentration")

#Calculate MAPE on Validation
test.hw_m=forecast(HWES.pm2_m,h=6)
error=pm2_valid-test.hw_m$mean
MAE=mean(abs(error))
MAPE=mean(abs(error)/abs(pm2_valid))
#MAE:2.029253
#MAPE: 0.2090463


#Calculate RMSE on Validation
error_sq=(pm2_valid-test.hw_m$mean)^2
MSE = mean(error_sq)
RMSE = sqrt(MSE)
#RMSE: 2.106654


