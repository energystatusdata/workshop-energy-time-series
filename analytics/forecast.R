
library(tseries)
library(forecast)

#### Create a list of train and test data sets for a particular feature
splitData <- function(d, feature){
  res <- list()
  for(i in 1:12){
    days <- unique(format(d[as.integer(format(d$timestamp, format="%m")) == i,]$timestamp, format="%d"))
    days <- days[days > 21]
    for(day in days) {
      dtrain <- d[(as.integer(format(d$timestamp, format="%m")) == i) & (as.integer(format(d$timestamp, format="%d")) < day), feature]
      dtest <- d[(as.integer(format(d$timestamp, format="%m")) == i) & (as.integer(format(d$timestamp, format="%d")) == day), feature]
      res[[length(res)+1]] <- list(dtrain, dtest)
    }
  }
  return(res)
}

################# Exponential Smoothing

getExpSmooth <- function(traintest){
  expsm1 <- expsm2 <- expsm3 <- expsm4 <- expsm5 <- matrix(nrow = 24, ncol = length(traintest))
  i <- 1
  for(d1 in traintest){
    dtrain <- d1[[1]]
    
    tstrain <- ts(dtrain, frequency = 24)
    
    # The first letter in the error type, the second is the trend type and the thrid the season type
    # N = None
    # A = Additive
    fit <- ets(tstrain, "ANN")
    expsm1[, i] <- forecast(fit, h = 24, method ='ets')$mean
    fit <- ets(tstrain, "AAN")
    expsm2[, i] <- forecast(fit, h = 24, method ='ets')$mean
    fit <- ets(tstrain, "ANA")
    expsm3[, i] <- forecast(fit, h = 24, method ='ets')$mean
    fit <- ets(tstrain, "AAA")
    expsm4[, i] <- forecast(fit, h = 24, method ='ets')$mean
    
    print(paste0("Processing test sample ", i))
    i <- i + 1
  }
  list(expsm1, expsm2, expsm3, expsm4, expsm5)
}

forecast_ExpSmooth <- function(d, var.name){
  traintest <- splitData(d, var.name)
  test <- t(do.call(rbind, sapply(traintest,function(x) x[2])))
  
  forecasts <- getExpSmooth(traintest)
  
  errors <-   getErrors(forecasts, test)
  result <- matrix(ncol = length(forecasts) , nrow = 2)
  for(i in 1:length(errors)){
    result[i,] <- errors[[i]]
  }
  result <- as.data.frame(result, row.names = c("MSE", "RMSE"))
  colnames(result) <- c("ExpSmooth ANN", "ExpSmooth AAN", "ExpSmooth ANA", "ExpSmooth AAA")
  list(forecasts, test, result)
}

################# 

################# ARIMA

getAutoARIMA <- function(traintest){
  arima <- matrix(nrow = 24, ncol = length(traintest))
  models <- rep(NA, length(traintest))
  i <- 1
  for(d1 in traintest){
    dtrain <- d1[[1]]
    tstrain <- ts(dtrain, frequency = 24)
    fit <- auto.arima(tstrain)
    arima[, i] <- forecast(fit, h = 24)$mean
    models[i] <- forecast(fit, h = 24)$method
    print(paste0("auto ARIMA", i))
    i <- i + 1
  }
  list(arima, models)
}

forecast_ARIMA <- function(d, var.name){
  traintest <- splitData(d, var.name)
  test <- t(do.call(rbind, sapply(traintest,function(x) x[2])))
  
  ARIMA_results <- getAutoARIMA(traintest)
  forecasts <- ARIMA_results[1]
  models <- ARIMA_results[[2]]
  
  errors <-   getErrors(forecasts, test)
  result <- matrix(ncol = length(forecasts) , nrow = 2)
  for(i in 1:length(errors)){
    result[i,] <- errors[[i]]
  }
  result <- as.data.frame(result, row.names = c("MSE", "RMSE"))
  colnames(result) <- c("ARIMA")
  list(forecasts, test, models, result)
}

#################

################# SOME NAIVE BASELINES

# Take the last day and return its average in a 24-value vector
naiveLastDayAverage <- function(dts){
  x <- as.vector(dts[(length(dts) - 23):length(dts)])
  lastday <- (ts(x, frequency = 24, start = (length(dts)/24) + 1))
  return(rep(mean(lastday), 24))
}

# 1) naive function (package "forecast") - the last observation (naive1)
# 2) sanive() function (package "forecast"), returns the preceeding day (naive2)
# 3) naiveLastDayAverage (preceeding day), returns the average of the preceeding day (naive3)
getNaive <- function(traintest){
  naive1 <- naive2 <- naive3 <- actual <- matrix(nrow = 24, ncol = length(traintest))
  
  i <- 1
  for(d1 in traintest){
    dtrain <- d1[[1]]
    dtest <- d1[[2]]
  
    tstrain <- ts(dtrain, frequency = 24)
    tstest <- ts(dtest, frequency = 24, start = (length(tstrain)/24) + 1)
    
    naive1[, i] <- naive(tstrain, h = 24)$mean
    naive2[, i] <- snaive(tstrain, h = 24)$mean
    naive3[, i] <- naiveLastDayAverage(tstrain)
    print(paste0("naive", i))
    i <- i + 1
  }
  list(naive1, naive2, naive3)
}

forecast_Naive <- function(d, var.name){
  traintest <- splitData(d, var.name)
  test <- t(do.call(rbind, sapply(traintest,function(x) x[2])))
  
  forecasts <- getNaive(traintest)
  
  errors <-   getErrors(forecasts, test)
  result <- matrix(ncol = length(forecasts) , nrow = 2)
  for(i in 1:length(errors)){
    result[i,] <- errors[[i]]
  }
  result <- as.data.frame(result, row.names = c("MSE", "RMSE"))
  colnames(result) <- c("naive1", "naive2", "naive3")
  list(forecasts, test, result)
}

#################

################# ERROR COMPUTATION AND PLOTTING

getErrors <- function(forecasts, test){
  MSE <- RMSE <- rep(NA, length(forecasts))
  for(i in 1:(length(forecasts))){
    MSE[i] <- mean((forecasts[[i]] - test)^2)
    RMSE[i] <- sqrt(MSE[[i]])
  }
  return(list(MSE, RMSE))
}

par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.7) # just reduce the padding
plot_forecast <- function(reference, forecast, title) {
  layout(matrix(c(1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13), 4, 4, byrow = TRUE))
  plot(as.vector(reference), type="l", ann=FALSE)
  title(main=title)
  lines(as.vector(forecast), col="red")
  
  names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aou", "Sep", "Oct", "Nov", "Dec")
  i <- 1
  for(x in c(9,18,28,37,47,56,66,76,85,95,104,113)) { # Plot each the day at the end of each month
    plot(as.vector(reference[(x*24):(x*24+24)]), type="l", xaxt='n', ann=FALSE, ylim=c(0, 30000))
    title(main = names[i])
    lines(as.vector(forecast[(x*24):(x*24+24)]), col="red")
    i <- i + 1
  }
}

################# 

###########################################################################

############################# OPEN AND PREPROCESS THE DATA
d <- read.csv("data/weather_and_power_generation.csv")

d[is.na(d[,"DE_wind_generation"]), "DE_wind_generation"] <- mean(d[,"DE_wind_generation"], na.rm = TRUE)
d[is.na(d[,"DE_solar_generation"]), "DE_solar_generation"] <- mean(d[,"DE_solar_generation"], na.rm = TRUE)

###################### Timestamp preprocessing
d$timestamp <- as.POSIXct(d$X, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
d$X <- NULL
######################

##############################
  
##################### Run Exponential Smoothing

wind_ExpSmooth <- forecast_ExpSmooth(d, "DE_wind_generation")
wind_ExpSmooth
#     ExpSmooth ANN ExpSmooth AAN ExpSmooth ANA ExpSmooth AAA
#MSE   18611941.304  17893873.898    18663267.2  18692672.507
#RMSE      4314.156      4230.115        4320.1      4323.502

plot_forecast(wind_ExpSmooth[[2]], wind_ExpSmooth[[1]][[1]], title=paste("Exponential Smoothing ANN - Wind (RMSE:", wind_ExpSmooth[[3]][2,1], ")", sep=""))
plot_forecast(wind_ExpSmooth[[2]], wind_ExpSmooth[[1]][[2]], title=paste("Exponential Smoothing AAN - Wind (RMSE:", wind_ExpSmooth[[3]][2,2], ")", sep=""))
plot_forecast(wind_ExpSmooth[[2]], wind_ExpSmooth[[1]][[3]], title=paste("Exponential Smoothing ANA - Wind (RMSE:", wind_ExpSmooth[[3]][2,3], ")", sep=""))
plot_forecast(wind_ExpSmooth[[2]], wind_ExpSmooth[[1]][[4]], title=paste("Exponential Smoothing AAA - Wind (RMSE:", wind_ExpSmooth[[3]][2,4], ")", sep=""))

solar_ExpSmooth <- forecast_ExpSmooth(d, "DE_solar_generation")
solar_ExpSmooth
#     ExpSmooth ANN ExpSmooth AAN ExpSmooth ANA ExpSmooth AAA
#MSE   52966748.961  52966749.052   3634951.447   3679579.134
#RMSE      7277.826      7277.826      1906.555      1918.223

plot_forecast(solar_ExpSmooth[[2]], solar_ExpSmooth[[1]][[1]], title=paste("Exponential Smoothing ANN - Solar (RMSE:", solar_ExpSmooth[[3]][2,1], ")", sep=""))
plot_forecast(solar_ExpSmooth[[2]], solar_ExpSmooth[[1]][[2]], title=paste("Exponential Smoothing AAN - Solar (RMSE:", solar_ExpSmooth[[3]][2,2], ")", sep=""))
plot_forecast(solar_ExpSmooth[[2]], solar_ExpSmooth[[1]][[3]], title=paste("Exponential Smoothing ANA - Solar (RMSE:", solar_ExpSmooth[[3]][2,3], ")", sep=""))
plot_forecast(solar_ExpSmooth[[2]], solar_ExpSmooth[[1]][[4]], title=paste("Exponential Smoothing AAA - Solar (RMSE:", solar_ExpSmooth[[3]][2,4], ")", sep=""))

##################### Run ARIMA

wind_ARIMA <- forecast_ARIMA(d, "DE_wind_generation")
wind_ARIMA
#            ARIMA
#MSE  18253698.252
#RMSE     4272.435

plot_forecast(wind_ARIMA[[2]], wind_ARIMA[[1]][[1]], title=paste("ARIMA - Wind (RMSE:", wind_ARIMA[[4]][2,1], ")", sep=""))
table(wind_ARIMA[[3]]) # just watch the frequency of ARIMA parameters

solar_ARIMA <- forecast_ARIMA(d, "DE_solar_generation")
solar_ARIMA
#           ARIMA
#MSE  8458544.409
#RMSE    2908.358

plot_forecast(solar_ARIMA[[2]], solar_ARIMA[[1]][[1]], title=paste("ARIMA - Solar (RMSE:", solar_ARIMA[[4]][2,1], ")", sep=""))
table(solar_ARIMA[[3]])# just watch the frequency of ARIMA parameters

##################### Run Baselines

wind_Naive <- forecast_Naive(d, "DE_wind_generation")
wind_Naive
#           naive1       naive2       naive3
#MSE  18611759.204 37058763.443 30213389.321
#RMSE     4314.135     6087.591     5496.671

plot_forecast(wind_Naive[[2]], wind_Naive[[1]][[1]], title=paste("Naive 1 - Wind (RMSE:", wind_Naive[[3]][2,1], ")", sep=""))
plot_forecast(wind_Naive[[2]], wind_Naive[[1]][[2]], title=paste("Naive 2 - Wind (RMSE:", wind_Naive[[3]][2,2], ")", sep=""))
plot_forecast(wind_Naive[[2]], wind_Naive[[1]][[3]], title=paste("Naive 3 - Wind (RMSE:", wind_Naive[[3]][2,3], ")", sep=""))

solar_Naive <- forecast_Naive(d, "DE_solar_generation")
solar_Naive
#           naive1      naive2       naive3
#MSE  52966748.961 4190024.051 32280177.310
#RMSE     7277.826    2046.955     5681.565

plot_forecast(solar_Naive[[2]], solar_Naive[[1]][[1]], title=paste("Naive 1 - Solar (RMSE:", solar_Naive[[3]][2,1], ")", sep=""))
plot_forecast(solar_Naive[[2]], solar_Naive[[1]][[2]], title=paste("Naive 2 - Solar (RMSE:", solar_Naive[[3]][2,2], ")", sep=""))
plot_forecast(solar_Naive[[2]], solar_Naive[[1]][[3]], title=paste("Naive 3 - Solar (RMSE:", solar_Naive[[3]][2,3], ")", sep=""))

##################### 
