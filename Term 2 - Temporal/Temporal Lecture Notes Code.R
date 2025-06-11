#### Lecture Notes R Code ####


## 1. Introduction to Time Series
# Using library 'astsa'
# install.packages('astsa')
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



### 2. Linear Systems

# Linear recurrence relation
library(astsa)
tsplot(0:50, sapply(0:50, function(t) 10*(0.9)^t),
       type="p", pch=19, cex=0.5, col=4, ylab="y_t",
       main="y_0 = 10, phi = 0.9")

# In continuous time
tsplot(0:50, sapply(0:50, function(t) 10*exp(-0.1*t)),
       cex=0.5, col=4, ylab="y(t)", main="y_0 = 10, lambda = 0.1")

# AR(1) model
tsplot(filter(rnorm(80, 0, 0.5), 0.9, "rec", init=10),
       type="p", col=4, pch=19, cex=0.5, ylab="Y_t",
       main="y_0=10, alpha=0.9, sigma=0.5")
abline(h=0, col=3, lwd=2)

# Continuous linear Markov process - Ornstein-Uhlenbeck process
T=50; dt=0.1; lambda=0.2; sigma=0.5; y0=10
times = seq(0, T, dt)
y = filter(rnorm(length(times), 0, sigma*sqrt(dt)), 
           1 - lambda*dt, "rec", init=y0)
tsplot(times, y, col=4, ylab="Y(t)", 
       main=paste("lambda =",lambda,"sigma =",sigma,"y0 =",y0))
abline(h=0, col=3, lwd=2)

# Second order linear recursion
tsplot(filter(rep(0,50), c(1.5, -0.75), "rec", init=c(10,10)),
       type="p", pch=19, col=4, cex=0.5, ylab="Y_t",
       main="phi_1 = 1.5, phi_2 = -0.75")
abline(h=0, col=3, lwd=2)

# AR(2) mode
tsplot(filter(rnorm(80), c(1.5, -0.75), "rec", init=c(20,20)),
       type="p", pch=19, col=4, cex=0.5, ylab="Y_t",
       main="phi_1 = 1.5, phi_2 = -0.75")
abline(h=0, col=3, lwd=2)



### 3. ARMA Models

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
polyroot(c(1, -1.5, 0.75))

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



### 4. Estimating and Forecasting

## 4.1 Fitting models

# AR(2) moment matching
tsplot(sunspot.year, col=4, main='Annual Sunspot Data')
tsplot(sqrt(sunspot.year), col=4)

# Detrend first
x <- detrend(sqrt(sunspot.year))
tsplot(x, col=4)

# Check periodicity by ACF
acf2(x)

# ACF suggests period of 10-11years; PACF suggests AR(2), which we fit
rho = acf(x, plot=FALSE, lag.max=2)$acf[2:3]
P = matrix(c(1, rho[1], rho[1], 1), ncol=2)
phi = solve(P, rho)
phi

# These coefficients correspond to complex roots of the char. polynomial - confirm period
u = polyroot(c(1, -phi))[1]
f = Arg(u)/(2*pi)
1/f

# AR(2) Least Squares
n = length(x)
X = cbind(x[2:(n-1)], x[1:(n-2)])
xp = x[3:n]
lm(xp ~ 0 + X)

# Fitting ARMA(2,1) model LS
loss = function(param) {
  phi = param[1:2]; theta = param[3]
  eps = signal::filter(c(1, -phi), c(1, theta),
                       x[3:length(x)], init.x=x[1:2], init.y=c(0))
  sum(eps*eps)
}

optim(rep(0.1,3), loss)$par

# Instead using arima function
arima(x, c(2, 0, 1), include.mean=FALSE)


## 4.2 Forecasting with models

# Forecasting AR(2)
n = length(x); k=50         # k=50 forecasts
phi = arima(x, c(2,0,0))$coef[1:2]       # fit the model and extract coefs
fore = filter(rep(0,k), phi, "rec", init=x[n:(n-1)])    # forecast by fitting the ar(2) with initial values n,n-1
tsplot(x, xlim=c(tsp(x)[1], tsp(x)[2]+k), col=4)      # tsp(x) gives (start, end, freq.)
lines(seq(tsp(x)[2]+1, tsp(x)[2]+k, frequency(x)), fore, 
      col=2, lwd=2)

# With variance
mod = arima(x, c(2,0,0))
fore = predict(mod, n.ahead=k)
pred = fore$pred
sds = fore$se
tsplot(x, xlim=c(tsp(x)[1], tsp(x)[2]+k),
       ylim=c(-10, 10), col=4)
ftimes = seq(tsp(x)[2]+1, tsp(x)[2]+k, tsp(x)[3])
lines(ftimes, pred, col=2, lwd=2)
lines(ftimes, pred+2*sds, col=2)
lines(ftimes, pred-2*sds, col=2)



### 5. Spectral Analysis
library(astsa)

nu <- seq(-0.5, 0.5, 0.001)

specd <- function(nu) {
  z <- complex(1, cos(2*pi*nu), -sin(2*pi*nu))
  phi <- 1 - 1.5*z + 0.75*z*z
  1/abs(phi)^2
}

spec <- sapply(nu, specd)

tsplot(ts(spec, start=-0.5, deltat=0.001),
       col=3, lwd=1.5, xlab='nu', ylab='spectral density')

# Since even, only plot half the function
nu = seq(0, 0.5, 0.002)
spec = sapply(nu, specd)
tsplot(ts(spec, start=0, deltat=0.002),
       col=3, lwd=1.5, xlab="nu", ylab="spectral density")

# On a half-log scale
tsplot(ts(log(spec), start=0, deltat=0.002),
       col=3, lwd=1.5, xlab='frequency (nu)', ylab='log(spectral density)')

# Finding the peak and oscillation period
peak = which.max(spec)*0.002
1/peak    # = 12 so oscillations with period 12

# Function to invert the DFT
ifft = function(x) fft(x, inverse=TRUE) / length(x)

# Checking the function works
x = 1:5
fx= fft(x)
ifft(fx)
Re(ifft(fx))

# Periodogram for AR(2) model
# First simulating some data
set.seed(42)
n=500
x = arima.sim(n=n, list(ar = c(1.5, -0.75)))
tsplot(x, col=4)

# To plot the basic periodogram
tsplot(abs(fft(x))^2, col=3, lwd=1.5, xlab='k', ylab='Power')

# Since symmetric, only plot the first half
tsplot(abs(fft(x)[2:(n/2+1)])^2, col=3, lwd=1.5, xlab='k', ylab='Power')

# On a half-log scale
tsplot(log(abs(fft(x)[2:(n/2+1)])^2), col=3, lwd=1.5, xlab='k', ylab='log(Power)')

# Like R's built periodogram function, which does a little pre-processing
spec.pgram(x, taper=0, detrend=FALSE, col=3, lwd=1.5, main='Raw periodogram')

spec.pgram(x, col=3, lwd=1.5, main='Raw periodogram')

spectrum(x, spans=c(5,7), col=3, lwd=2, main='Smoothed periodogram')

# Monthly sunspot dataset
tsplot(sunspot.month, col=4)

# Periodogram
spectrum(sunspot.month, spans=c(5,7), col=3)
# Shows strong peak at low frequency, but a lot of higher frequency behaviour too

# Checking the DFT
ft = fft(sunspot.month)
tsplot(log(abs(ft)), col=3, xlab='Frequency', ylab='log(DFT)')
# Interested in peaks at end, want to remove noise in the middle - keep first and last
# 150 frequency components
ft[152:(length(ft)-152+2)]=0
tsplot(log(abs(ft)), col=3, xlab='Frequency', ylab='log(DFT)')

# Switching back to time domain
sm = Re(ifft(ft))
tsplot(sm, col=4, lwd=1.5, ylab='Smoothed monthly sunspot data')



### 6. Hidden Markov Models

library(astsa)
y = BCJ[,"boa"]
length(y)

tsplot(y, col=4, main="BoA returns")

# HMM to segment ts into regions of high and low volatility
hist(y[2000:3000], 30, col=4, freq=FALSE,
     main="Region of low volatility")
curve(dnorm(x, 0, 0.015), add=TRUE, col=3, lwd=2)

hist(y[1000:1500], 30, col=4, freq=FALSE,
     main="Region of high volatility")
curve(dcauchy(x, 0, 0.025), add=TRUE, col=3, lwd=2)

# Implementing a filter
hmmFilter = function(P, l)
  function(f, y) {
    fNew = (f %*% P) * l(y)
    fNew / sum(fNew)
  }

advance = hmmFilter(
  matrix(c(0.999, 0.001, 0.005, 0.995), ncol=2, byrow=TRUE),
  function(y) c(dnorm(y, 0, 0.015), dcauchy(y, 0, 0.025))
)

Reduce(advance, y, c(0.5, 0.5))

# Plot full set of filtered probabilities
fpList = Reduce(advance, y, c(0.5, 0.5), acc=TRUE)
fpMat = sapply(fpList, cbind)
fp2Ts = ts(fpMat[2, -1], start=start(y), freq=frequency(y))
tsplot(y, col=4, ylim=c(-0.3, 1))
lines(fp2Ts, col=2, lwd=1.5)

# Modify function to update marginal likelihood of data so far, in addition to filtered probabilities
hmmFilterML = function(P, l)
  function(fl, y) {
    fNew = (fl$f %*% P) * l(y)
    ml = sum(fNew)
    list(f=fNew/ml, ll=(fl$ll + log(ml)))
  }

advance = hmmFilterML(
  matrix(c(0.999, 0.001, 0.005, 0.995), ncol=2, byrow=TRUE),
  function(y) c(dnorm(y, 0, 0.015), dcauchy(y, 0, 0.025))
)

Reduce(advance, y, list(f=c(0.5, 0.5), ll=0))

# Now smoothing
hmmSmoother = function(P)
  function(fp, sp) {
    fp * ((sp / (fp %*% P)) %*% t(P))
  }

backStep = hmmSmoother(
  matrix(c(0.999, 0.001, 0.005, 0.995), ncol=2, byrow=TRUE)
)

spList = Reduce(backStep, fpList, right=TRUE, acc=TRUE)
spMat = sapply(spList, cbind)
sp2Ts = ts(spMat[2, -1], start=start(y), freq=frequency(y))
tsplot(y, col=4, ylim=c(-0.3, 1))
lines(sp2Ts, col=2, lwd=1.5)

# Generating samples
hmmSampler = function(P) {
  p = nrow(P)
  function(fp, x) {
    sample(1:p, 1, prob=fp*P[,x])
  }
}

backSample = hmmSampler(
  matrix(c(0.999, 0.001, 0.005, 0.995), ncol=2, byrow=TRUE)
)

set.seed(42)
xList = Reduce(backSample, head(fpList, -1),
               init=sample(1:2, 1, prob=tail(fpList, 1)[[1]]),
               right=TRUE, acc=TRUE)
xVec = unlist(xList)
xTs = ts(xVec[-1], start=start(y), freq=frequency(y))
tsplot(y, col=4, ylim=c(-0.3, 1))
lines(xTs-1, col=2, lwd=1.5)

# Parameter estimation
mll = function(theta) {
  advance = hmmFilterML(
    matrix(c(0.999, 0.001, 0.005, 0.995), ncol=2, byrow=TRUE),
    function(yt) 
      c(dnorm(yt, 0, theta[1]), dcauchy(yt, 0, theta[2]))
  )
  Reduce(advance, y, list(f=c(0.5, 0.5), ll=0))$ll
}

mll(c(0.015, 0.025))

optim(c(0.015, 0.025), mll, control=list(fnscale=-1))



### 7. Dynamic Linear Models

# Filtering
kFilter = function(G, W, F, V)
  function(mC, y) {
    m = G %*% mC$m
    C = (G %*% mC$C %*% t(G)) + W
    f = F %*% m
    Q = (F %*% C %*% t(F)) + V
    K = t(solve(Q, F %*% C))
    m = m + (K %*% (y - f))
    C = C - (K %*% Q %*% t(K))
    list(m=m, C=C)
  }

library(astsa)
tsplot(soi, col=4)

# Suppose SRW appropriate so G=1, and observations are noisy observations of hidden trend so F=1.
# Also suppose monthyl change in long run trend has sd=0.01, and noise variance for observations has sd=0.5
advance = kFilter(1, 0.01^2, 1, 0.5^2)
Reduce(advance, soi, list(m=0, C=100))

# Keeping full set of filtered states and plot over data
fs = Reduce(advance, soi, list(m=0, C=100), acc=TRUE)
fsm = sapply(fs, function(s) s$m)
fsTs = ts(fsm[-1], start=start(soi), freq=frequency(soi))
tsplot(soi, col=4)
lines(fsTs, col=2, lwd=2)

# Modify filter to compute marginal likelihood
kFilterML = function(G, W, F, V)
  function(mCl, y) {
    m = G %*% mCl$m
    C = (G %*% mCl$C %*% t(G)) + W
    f = F %*% m
    Q = (F %*% C %*% t(F)) + V
    ll = mvtnorm::dmvnorm(y, f, Q, log=TRUE)
    K = t(solve(Q, F %*% C))
    m = m + (K %*% (y - f))
    C = C - (K %*% Q %*% t(K))
    list(m=m, C=C, ll=mCl$ll+ll)
  }

advance = kFilterML(1, 0.01^2, 1, 0.5^2)
Reduce(advance, soi, list(m=0, C=100, ll=0))

# Parameter estimation
# Estimating variance param. V and W
mll = function(wv) {
  advance = kFilterML(1, wv[1], 1, wv[2])
  Reduce(advance, soi, list(m=0, C=100, ll=0))$ll
}

mll(c(0.01^2, 0.5^2))

optim(c(0.01^2, 0.5^2), mll, control=list(fnscale=-1))

# Using the 'dlm' package
# Fitting the DLM to data
library(dlm)
mod = dlm(FF=1, GG=1, V=0.5^2, W=0.01^2, m0=0, C0=100)
fs = dlmFilter(soi, mod)
tsplot(soi, col=4, main="Filtered states")
lines(fs$m, col=2, lwd=2)

# Smoothed states
ss = dlmSmooth(soi, mod)
tsplot(soi, col=4, main="Smoothed states")
lines(ss$s, col=2, lwd=2)

# Optimise param
buildMod = function(lwv)
  dlm(FF=1, GG=1, V=exp(lwv[2]), W=exp(lwv[1]),
      m0=0, C0=100)

opt = dlmMLE(soi, parm=c(log(0.5^2), log(0.01^2)),
             build=buildMod)
opt
exp(opt$par)

# Smoothed states for optimised model
mod = buildMod(opt$par)
ss = dlmSmooth(soi, mod)
tsplot(soi, col=4, main="Smoothed states")
lines(ss$s, col=2, lwd=1.5)



### 8. State Space Modelling

library(astsa)
library(dlm)
y = log(jj)
tsplot(y, col=4, lwd=2, main="Log of J&J earnings")

# Creating model
mod = dlmModPoly(2, dV=0.1^2, dW=c(0.01^2, 0.01^2))       # 2 for locally linear, diagonal of V,W
mod

ss = dlmSmooth(y, mod)
tsplot(y, col=4, lwd=2, main="Log of J&J earnings")
lines(ss$s[,1], col=2, lwd=2)

# Fit model with locally linear trend, and quarterly seasonal effect
mod = dlmModPoly(2, dV=0.1^2, dW=c(0.01^2, 0.01^2)) +
  dlmModSeas(4, dV=0, dW=c(0.02^2, 0, 0))
mod

ss = dlmSmooth(y, mod) # smoothed states
so = ss$s[-1,] %*% t(mod$FF) # smoothed observations
so = ts(so, start=start(y), freq=frequency(y))
tsplot(y, col=4, lwd=2, main="Log of J&J earnings")
lines(ss$s[,1], col=3, lwd=1.5)
lines(so, col=2, lwd=1.2)

# Now for forecasting
fit = dlmFilter(y, mod)
fore = dlmForecast(fit, 16) # forecast 16 time points
pred = ts(c(tail(y, 1), fore$f), 
          start=end(y), frequency=frequency(y))
upper = ts(c(tail(y, 1), fore$f + 2*sqrt(unlist(fore$Q))), 
           start=end(y), frequency=frequency(y))
lower = ts(c(tail(y, 1), fore$f - 2*sqrt(unlist(fore$Q))), 
           start=end(y), frequency=frequency(y))
all = ts(c(y, upper[-1]), start=start(y),
         frequency = frequency(y))
tsplot(all, ylab="log(earnings)",
       main="Forecasts with 2SD intervals")
lines(y, col=4, lwd=1.5)
lines(pred, col=2, lwd=2)
lines(upper, col=2)
lines(lower, col=2)

# Now monthly birth data
y = birth
tsplot(y, col=4, lwd=1.5, main="Monthly births")

# Seasonal model on 2 Fourier harmonics, estimate varianecs with ML
buildMod = function(lpar)
  dlmModPoly(1, dV=exp(lpar[1]), dW=exp(lpar[2])) +
  dlmModTrig(12, 2, dV=0, dW=rep(exp(lpar[3]), 4))
opt = dlmMLE(y, parm=log(c(100, 1, 1)),
             build=buildMod)
opt

mod = buildMod(opt$par)
ss = dlmSmooth(y, mod)
so = ss$s[-1,] %*% t(mod$FF) # smoothed observations
so = ts(so, start(y), frequency = frequency(y))
tsplot(y, col=4, lwd=1.5, main="Monthly births")
lines(ss$s[,1], col=3, lwd=1.5)
lines(so, col=2, lwd=1.2)

# Locally constant model, with seasonal effect, and residual correlation structure with an AR(2)
y = soi
buildMod = function(param)
  dlmModPoly(1, dV=exp(param[1]), dW=exp(param[2])) +
  dlmModARMA(ar=c(param[3], param[4]), 
             sigma2=exp(param[5]), dV=0) +
  dlmModTrig(12, 2, dV=0, dW=rep(exp(param[6]), 4))
opt = dlmMLE(y, parm=c(log(0.1^2), log(0.01^2), 
                       0.2, 0.1, log(0.1^2), log(0.01^2)), build=buildMod)
opt

mod = buildMod(opt$par)
ss = dlmSmooth(y, mod)
so = ss$s[-1,] %*% t(mod$FF) # smoothed observations
so = ts(so, start(y), freq=frequency(y))
tsplot(y, col=4, lwd=1.5, main="SOI")
lines(ss$s[,1], col=3, lwd=1.5)
lines(so, col=2)


