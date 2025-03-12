### Introduction to Time Series ###

## Question 1:

# Using library 'astsa'
install.packages('astsa')
library(astsa)

# Chicken dataset
plot(chicken)

# Time series plot of chicken data
tsplot(chicken)
tsplot(chicken, col=4, lwd=2, 
       main="Price of chicken", ylab="US cents per pound")

# Detrending by removing linear fit
tc = time(chicken)
res = residuals(lm(chicken ~ tc))
dc = ts(res, start=start(chicken), freq=frequency(chicken))
plot(dc, lwd=2, col=4, main="Detrended price of chicken")

# Using detrend function
tsplot(detrend(chicken), lwd=2, col=4)

# Autocorrelation function
lag1.plot(dc, pch=19, col=4)
acf1(dc)

# Moving average filter
tsplot(chicken, col=4,
       main="Chicken price smoothed with a moving average filter")
lines(filter(chicken, rep(1/25, 25)), col=2, lwd=2)

# Quarterly share earnings for J&J
tsplot(jj, col=4, lwd=2, main="Quarterly J&J earnings")

# Plot of the log for linear trend
tsplot(log(jj), col=4, lwd=2, main="Log of J&J earnings")

# Smoothed with linear filter
sjj = filter(log(jj), c(0.125,0.25,0.25,0.25,0.125))
tsplot(sjj, col=4, lwd=2, main="Smoothed log of J&J earnings")

# Detrend by removing the smoothed linear filter
dljj = log(jj) - sjj
tsplot(dljj, col=4, lwd=2, main="Detrended log of J&J earnings")

# Directly computing detrended data
dljj = filter(log(jj), c(-0.125,-0.25,0.75,-0.25,-0.125))
tsplot(dljj, col=4, lwd=2, main="Detrended log of J&J earnings")

# Estimate of seasonal effects of each quarter
seff = rowMeans(matrix(dljj, nrow=4), na.rm=TRUE)
seff

# Stripping seasonal effects
dsdtljj = dljj - seff
tsplot(dsdtljj, col=4, lwd=2, 
       main="Detrended and deseasonalised log of J&J earnings")

# Using library 'signal'
library(signal)

# Exponential smoothing
alpha = 0.1
sm = signal::filter(alpha, c(1, alpha-1), chicken, init.y=chicken[1])
tsp(sm) = tsp(chicken)
tsplot(chicken, col=4)
lines(sm, col=2, lwd=2)

# Multivariate time series
str(climhyd)
head(climhyd)

tsplot(climhyd, ncol=2, col=2:7)

pairs(climhyd, col=4, pch=19, cex=0.5)



## Question 2

# Using cardox dataset showing carbon dioxide levels at a location in Hawaii
cardox
tsplot(cardox, lwd=2, col=4, main='Mean Carbon Dioxide Level')

# Detrend by applying convolutional linear filter
f = c(1/24, rep(1/12, 11), 1/24)
scardox = filter(cardox, f)

tsplot(scardox, col=4, lwd=2, main = 'Deaseasonalised Mean Carbon Dioxide Level')

# Estimating monthly seasonal effects
dtcardox <- cardox- scardox
se <- rowMeans(matrix(dtcardox, nrow=12), na.rm=TRUE)

barplot(se, names.arg=c(month.name[3:12], month.name[1:2]), col=4)

## Question 3

# Second order Markovian scalar linear Gaussian auto-regressive model - AR(2)
set.seed(42)
tseries <- filter(rnorm(1000), c(1.5, -0.75), "rec", init=c(0,0))
tsplot(tseries,
       type="p", pch=19, col=4, cex=0.5, ylab="Y_t",
       main="phi_1 = 1.5, phi_2 = -0.75")
abline(h=0, col=3, lwd=2)

# Empirical mean and variance
mean(tseries)
var(tseries)

# Computing first few auto-correlations of the process
acf(tseries, lag.max=5, plot=FALSE)

# Reversing the time series and calculating empirical mean and variance
revtseries <- replace(tseries, TRUE, rev(tseries))

mean(revtseries)
var(revtseries)

# Computing first few auto-correlations of the process
acf(revtseries, lag.max=5, plot=FALSE)



## Question 4

# Simulating 10,000 observations from a centred VAR(1) process
phi <- matrix(c(0.9, 0.2, -0.1, 0.8), ncol=2, byrow=TRUE)
xList <- Reduce(function(xi, z) phi %*% xi + rnorm(2), 
               rep(0, 9999), c(0, 0), acc=TRUE)
x <- t(sapply(xList, cbind))
dim(x)

# Computing eigenvalues of Phi to confirm stability
abs(eigen(phi)$values)

# Computing empirical covariance matrix, and compare to its theoretical value
var(x)

library(netcontrol)
V <- dlyap(t(phi), diag(2))

# Compute a few empirical cross-correlations between the 2 components of the simulated time series
ccf(x[,1], x[,2], lag.max=5, plot=FALSE)

# Compare to theoretical values as suggested by stationary covariance function
acfList = Reduce(function(m, z) m %*% t(phi), 1:5, V, acc=TRUE) # Auto-covariance matrices
pos = sapply(acfList, function(m) m[2,1]) # cross-covariances for non-neg t
neg = sapply(acfList, function(m) m[1,2]) # for non-positive t
true_ccvf = c(rev(neg[-1]), pos) # cross-covariance function
true_ccf = true_ccvf / sqrt(V[1,1]*V[2,2]) # cross-correlation function
true_ccf

# Reverse time series and confirm empirical variance matrix unchanged
rx <- x[10000:1,]
var(rx)

# Investigate cross-correlations between components of reversed series
ccf(rx[,1], rx[,2], lag.max=5, plot=FALSE)











     