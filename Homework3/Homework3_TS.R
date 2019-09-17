#------------------------------------#
#            Homework 3              #
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
library(dyn)
#Ask Dr Simmons
library(urca)

# Saving File Locations and Uploading CSV File #
air_quality <- read.csv("/PM_2_5_Raleigh2.csv")

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

#Create Data Split---------------------------------------------------------------------------------------
pm2_train= subset(monthly_pm2, end=length(monthly_pm2)-6)
pm2_valid=subset(monthly_pm2,start=length(monthly_pm2)-5)

#Check for Stationarity: Random Walk---------------------------------------------------------------------------
#Augmented Dickey Fuller Test for Lag 0 - 2
ADF.Pvalues <- rep(NA, 3)
for(i in 0:2){
  ADF.Pvalues[i+1] <- adf.test(monthly_pm2, alternative = "stationary", k = i)$p.value
}

#Another way to do this using the 'urca' package
df <- ur.df(monthly_pm2, type = c("drift"), lags=2)
summary(df)

#Correction for  Random Walk (Assuming alpha = 0.01)-----------------------------------------------------------------------
#Difference the Data 
stationary_pm2 <- diff(monthly_pm2, differences= 1)
plot(stationary_pm2)

adf.test(monthly_pm2, alternative = 'stationary',k=2)
rw.drift=Arima(monthly_pm2,order=c(0,1,0))
summary(rw.drift)

#Correction for  Deterministic Trend (Assuming alpha = 0.05)----------------------------------------------------------------
#linear model of TS data with linear indicies 
dtrend <- lm(monthly_pm2 ~ c(1:length(monthly_pm2)))
plot(resid(dtrend))

#check if this is right
arima.dtrend=Arima(monthly_pm2,xreg=c(1:length(monthly_pm2)),order=c(0,0,0))

#plot residuals 
res.dtrend <- dtrend$residuals
res.a.dtrend <- arima.dtrend$residuals
plot(res.dtrend, res.a.dtrend)  


#Check for White Noise-----------------------------------------------------------------------------------
White.LB <- rep(NA, 10)
for(i in 1:10){
  White.LB[i] <- Box.test(arima.dtrend$residuals, lag = i, type = "Lj", fitdf = 1)$p.value
}

White.LB <- pmin(White.LB, 0.2)
barplot(White.LB, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags", ylim = c(0, 0.2))
abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")
