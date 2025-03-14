#### ARMA Models #####


## Question 1

# AR(p) Models
# AR(2) example from previous
library(astsa)
set.seed(42)
phi = c(1.5, -0.75); sigma = 1; n = 2000; burn = 100
x = filter(rnorm(n+burn, 0, sigma), phi, "rec")[-(1:burn)]
tsplot(x, col=4, cex=0.5, ylab="X_t",
       main=paste0("phi_1=", phi[1], ", phi_2=", phi[2]))
abline(h=0, col=3)

# Checking if eigenvalues of the generator of equivalent VAR(1) model have modulus less than 1
Phi = matrix(c(1.5,-0.75,1,0), ncol=2, byrow=TRUE)
evals = eigen(Phi)$values
abs(evals)

# Empirical ACF of simulated data
acf1(x)

# Computing true covariances by using VAR(1) representation and solve for stationary variance matrix
Sig = matrix(c(1,0,0,0), ncol=2)
V = netcontrol::dlyap(t(Phi), Sig)
V
# First row gives first 2 auto-covariances

# Recursively generate more
initAC = V[1,]
acvf = filter(rep(0,49), c(1.5,-0.75), 'rec', init=initAC)
acvf=c(initAC[1], acvf)
acvf[1:20]

# Turn them into auto-correlations and overlay on the empirical ACF of simulated dataset
acrf = acvf[-1]/acvf[1]
acrf[1:20]

acf1(x)
points(1:length(acrf), acrf, col=4, pch=19, cex=0.8)

# To compute just the ACF
ARMAacf(ar=c(1.5,-0.75), lag.max=10)

# Can simulated data instead using function arima.sim, instead of recursion of white noise
y = arima.sim(n=(n+burn), list(ar=phi), sd=sigma)[-(1:burn)]
tsplot(y, col=4)

# Now also computing the PACF
ARMAacf(ar=c(1.5, -0.75), lag.max=4)
ARMAacf(ar=c(1.5, -0.75), lag.max=4, pacf=TRUE)

# Seeing empirical ACF and PACF
acf2(x, max.lag=10)


# MA(q) Models
# Computing ACF and PACF of an MA(1) process with theta=0.8
ARMAacf(ma=c(0.8), lag.max=6)
ARMAacf(ma=c(0.8), lag.max=6, pacf=TRUE)

# Simulate data either by applying a convolutional filter to noise, or by arima.sim
x = filter(rnorm(2000), c(1,0.8))
x = arima.sim(n=2000, list(ma=c(0.8)), sd=1)
tsplot(x, col=4)

acf2(x, max.lag=6)


# ARMA(p,q) Models
# Consider an ARMA(1,1) with param. 0.8,0.6 and compute ACF and PACF
ARMAacf(ar=c(0.8), ma=c(0.6), lag.max=6)
ARMAacf(ar=c(0.8), ma=c(0.6), lag.max=6, pacf=TRUE)

# Simulating from model and looking at empirical statistics
x = signal::filter(c(1, 0.6), c(1, -0.8), rnorm(2000))
x = arima.sim(n=2000, list(ar=c(0.8), ma=c(0.6)), sd=1)
## either of the above approaches is fine
tsplot(x, col=4)

acf2(x, max.lag=6)




## Question 2
# ARMA(1,1) model parameters 0.8,0.6 - calculating ACF and PACF
library(astsa)
set.seed(42)
x <- signal::filter(c(1, 0.6), c(1,-0.8), rnorm(5000))
acf2(x, max.lag=8)

x <- arima.sim(n=5000,list(ar=c(0.8),ma=c(0.6)),sd=1)
acf2(x, max.lag=8)

# Fitting an AR(3) model using moment matching, then least squares, then built-in arima function
# Moment matching
acfs <- acf(x, plot=FALSE, lag.max=3)$acf
rho <- acfs[2:4]
P <- toeplitz(acfs[1:3])
phi <- solve(P, rho)
phi

# Least squares
n <- length(x)
X <- cbind(x[3:(n-1)], x[2:(n-2)], x[1:(n-3)])
xp <- x[4:n]
lm(xp ~ 0 + X)

# Arima
ar3 <- arima(x, c(3,0,0), include.mean = FALSE) 
ar3
ar3$coef

# Plot ACF and PACF of model fitted by Arima
# Fitted model
barplot(ARMAacf(ar=ar3$coef, lag.max=8))
barplot(ARMAacf(ar=ar3$coef, lag.max=8, pacf=TRUE))

# True model
barplot(ARMAacf(ar=c(0.8), ma=c(0.6), lag.max=8))
barplot(ARMAacf(ar=c(0.8), ma=c(0.6), lag.max=8, pacf=TRUE))

# Fit alternative models AR(2), AR(4), ARMA(1,1), ARMA(2,1)
ar2 = arima(x, c(2, 0, 0), include.mean=FALSE)
ar4 = arima(x, c(4, 0, 0), include.mean=FALSE)
arma11 = arima(x, c(1, 0, 1), include.mean=FALSE)
arma21 = arima(x, c(2, 0, 1), include.mean=FALSE)
