# Load libraries
rm(list=ls())
cat("\014")

library(tseries)
library(forecast)
library(reshape2)

china <- read.csv("china.csv")

imp <- ts(data = china$ch.imp, start = c(1984,1), frequency = 12)
exp <- ts(data = china$ch.exp, start = c(1984,1), frequency = 12)

par(mfrow=c(2,1))
plot(imp)
plot(exp)

exp_train <- window(exp, end = c(2002,12))
exp_test_orig <- window(exp, start = c(2003,1))

imp_train <- window(imp, end = c(2002,12))
imp_test_orig <- window(imp, start = c(2003,1))

exp_test <- melt(exp_test_orig)
imp_test <- melt(imp_test_orig)

par(mfrow = c(1,1))
plot(exp_train, main = "Chinese Export Figures (in 100 million USD)", ylab = "Exports", xlab = "Time")

hw.exp <- HoltWinters(x = exp_train, seasonal = "add", alpha = 0.2, beta = 0.7, gamma = 0.4) 

#TABULAR
predicted <- forecast(hw.exp,h = 72,level=95)
predicted <- data.frame(predicted)
predicted <- predicted$Point.Forecast
results <- cbind(predicted, exp_test)
RMSE <- sqrt(mean((results$predicted-results$value)^2))
RMSE
#GRAPHICAL
par(mfrow = c(2,1))
plot(hw.exp)
legend("topleft", legend = c("Observed", "Fitted"), lty = 1, col = c("black", "red"), cex = 0.5)
plot(forecast(hw.exp, h = 72, level = 95), ylab = "Exports")
legend("topleft", legend = c("Observed", "Predicted", "95% PI"), lty = 1, col = c("black", "blue", "grey"), cex = 0.5)

# SARIMA 

par(mfrow=c(1,1))
plot(exp_train)
par(mfrow=c(2,1))
plot(log(exp_train))
acf(log(exp_train), lag.max = 48)

# The raw time series is clearly not stationary. Try differencing once:
ndiffs(log(exp_train))
exp_train_d1 <- diff(log(exp_train))
plot(exp_train_d1, ylab = "log(export)")
acf(exp_train_d1, lag.max = 48)
adf.test(exp_train_d1)

# There still seems to be monthly seasonality (period = 12). Let's try differencing for that.
exp_train.12 <- diff(exp_train_d1, lag = 12)
plot(exp_train.12)
acf(exp_train.12, lag.max = 48)

# This looks good, so choose d=1, D=1, s=12. Let's check ACF and PACF to choose p, q, P, Q.
acf(exp_train.12, lag.max = 48)
pacf(exp_train.12, lag.max = 48)

# q <= 1, p <= 2, and Q <= 1, P <= 2. d = D = 1
m1 <- arima(log(exp_train), order = c(1,1,1), seasonal = list(order = c(2,1,1), period = 12), method = "CSS")
m2 <- arima(log(exp_train), order = c(1,1,1), seasonal = list(order = c(1,1,1), period = 12), method = "CSS")
m3 <- arima(log(exp_train), order = c(2,1,1), seasonal = list(order = c(2,1,1), period = 12), method = "CSS")
m4 <- arima(log(exp_train), order = c(2,1,1), seasonal = list(order = c(1,1,1), period = 12), method = "CSS")
m5 <- arima(log(exp_train), order = c(1,1,1), seasonal = list(order = c(0,1,0), period = 12), method = "CSS")

# Let's see which model auto.arima chooses as being optimal
auto.arima(log(exp_train), allowdrift = F)

m6 <- arima(log(exp_train), order = c(0,1,1), seasonal = list(order = c(2,0,0), period = 12), method = "CSS-ML")

# Goodness-of-fit metrics
RMSE1 <- sqrt(mean(m1$residuals^2))
RMSE2 <- sqrt(mean(m2$residuals^2))
RMSE3 <- sqrt(mean(m3$residuals^2))
RMSE4 <- sqrt(mean(m4$residuals^2))
RMSE5 <- sqrt(mean(m5$residuals^2))
RMSE6 <- sqrt(mean(m6$residuals^2))

sigma2 <- c(m1$sigma2,m2$sigma2,m3$sigma2,m4$sigma2,m5$sigma2,m6$sigma2)
loglik <- c(m1$loglik,m2$loglik,m3$loglik,m4$loglik,m5$loglik,m6$loglik)
AIC<- c(m1$aic,m2$aic,m3$aic,m4$aic,m5$aic,m6$aic)
RMSE <- c(RMSE1,RMSE2,RMSE3,RMSE4,RMSE5,RMSE6) 
d <- data.frame(sigma2,loglik,AIC,RMSE)
d[order(d$RMSE),]
# m4 (2,1,1 x 1,1,1) and m3 (2,1,1 x 2,1,1) seem to be optimal. Checking test RMSE for these models to get best model 

predicted_m4 <- forecast(m4, h = 72, level = 95)
predicted_m4 <- data.frame(predicted_m4)
predicted_m4 <- predicted_m4$Point.Forecast
predicted_m4 <- exp(predicted_m4)
results_m4 <- cbind(predicted_m4, exp_test)
results_m4
RMSE_m4 <- sqrt(mean((results_m4$predicted-results_m4$value)^2))
RMSE_m4

predicted_m3 <- forecast(m3, h = 72, level = 95)
predicted_m3 <- data.frame(predicted_m3)
predicted_m3 <- predicted_m3$Point.Forecast
predicted_m3 <- exp(predicted_m3)
results_m3 <- cbind(predicted_m3, exp_test)
RMSE_m3 <- sqrt(mean((results_m3$predicted-results_m3$value)^2))
RMSE_m3

predicted_m5 <- forecast(m5, h = 72, level = 95)
predicted_m5 <- data.frame(predicted_m5)
predicted_m5 <- predicted_m5$Point.Forecast
predicted_m5 <- exp(predicted_m5)
results_m5 <- cbind(predicted_m5, exp_test)
RMSE_m5 <- sqrt(mean((results_m5$predicted-results_m5$value)^2))
RMSE_m5

# Since test RMSE is lower for m5 (1x1x1, 0x1x0) we choose that as the optimal model

f<-forecast(m5, h=72, level=0.95)
l<-ts(f$lower, start = c(2003, 1), frequency = 12)  #95% PI LL
h<-ts(f$upper, start = c(2003, 1), frequency = 12) #95% PI UL
pred<-f$mean #predictions
par(mfrow=c(1,1))
plot(exp_train, xlim=c(1984,2012), ylim=c(0, 1500), main = "Chinese Exports (in Million USD)", ylab = "Exports", xlab = "Time")
abline(v = 2003, lwd = 2, col = "black")
points(exp(pred), type = "l", col = "red")
points(exp(l), type = "l", col = "grey")
points(exp(h), type = "l", col = "grey")
points(exp(f$fitted),type="l", col = "blue")
points(exp_test_orig,type="l", col = "green")
exp(pred)
exp_test_orig
legend("topleft", legend = c("Observed", "Fitted", "Predicted", "95% PI"), lty = 1, col = c("black", "green", "blue", "red"), cex = 0.5)

# We should check model diagnostics as we normally would with an ARMA model. Specifically, checking the heteroscedasticity 
# assumption can lead to an optimal choice of variance-stabilizing transformation

tsdiag(m5)

# VAR
library(vars)
VARselect(y = data.frame(exp_train, imp_train))

# p=3 appears appropriate. AIC reduces further but it is not that significant
# For the sake of parsimony let's choose p=3
m.var <- VAR(y = data.frame(exp_train,imp_train),p=6)
plot(m.var)

# Forecasting
predicted <- predict(m.var, n.ahead = 72, ci = 0.95)
predicted_var <- predicted$fcst$exp_train
rmse_var <- sqrt(mean(((data.frame(predicted_var)$fcst)-exp_test)^2))
rmse_var

