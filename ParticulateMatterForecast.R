#------------------------------------#
#            Homework 2              #
#          PM2.5 Forecast            #
#                                    #
#         Mehak Shamdasani           #
#------------------------------------#
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
file.dir <- "/Users/mehak/Desktop/MSA/FALL2020/TimeSeries/Homework2/"
air_quality <- read.csv("/Users/mehak/Desktop/MSA/FALL2020/TimeSeries/Homework2/PM_2_5_Raleigh2.csv")

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

#Create Data Split---------------------------------------------------------------------------


#Holt-Winters Exponential Smoothing Model---------------------------------------------------------------------------
#Additive
HWES.pm2 <- hw(monthly_pm2, seasonal = "additive")
summary(HWES.pm2)

plot(HWES.pm2, main = "Monthly PM2.5 Concentration with Holt-Winters ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 2008.25, col = "red", lty = "dashed")

autoplot(HWES.pm2)+
  autolayer(fitted(HWES.pm2),series="Fitted")+ylab("PM2.5 Concentration")

#Multiplicative
HWES.pm2 <- hw(monthly_pm2, seasonal = "multiplicative")
summary(HWES.pm2)

plot(HWES.pm2, main = "Monthly PM2.5 Concentration with Holt-Winters ESM Forecast", xlab = "Date", ylab = "PM2.5 Concentration")
#abline(v = 2008.25, col = "red", lty = "dashed")

autoplot(HWES.pm2)+
  autolayer(fitted(HWES.USAir),series="Fitted")+ylab("PM2.5 Concentration")

