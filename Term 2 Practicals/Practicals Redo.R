# Carbon dioxide time series needs deseasonalising and estimate seasonal effects
f = c(1/24, rep(1/12, 11), 1/24)
scardox = filter(cardox, f)
tsplot(scardox, lwd=2, col=4)

dtcardox = cardox-scardox
se = rowMeans(matrix(dtcardox, nrow=12), na.rm=TRUE)
se

# Simulate 1000 observations from AR(2) with phi = 1.5, -0.75, sigma=1
set.seed(42)
x = filter(rnorm(1000), c(1.5,-0.75), "rec", init=c(0, 0))

# Compute empirical mean and variance
mean(x)
var(x)

# Compute the first few auto-correlations
acf(x, lag.max=5, plot=FALSE)

# Simulate 5000 observations from ARMA(1,1) and check ACF and PACF
x=arima.sim(n=5000, list(ar=c(0.8), ma=c(0.6)), sd=1)
acf2(x, max.lag=8)

# Fit an AR(3) model by moment-matching, least squares and then arima
rho = acf(x, plot=FALSE, lag.max=3)$acf[2:4]
P = matrix(c(1, rho[1], rho[2],
           rho[1], 1, rho[1],
           rho[2], rho[1], 1), ncol=3)
phi = solve(P, rho)

n = length(x)
X = cbind(x[3:(n-1)], x[2:(n-2)], x[1:(n-3)])
xp = x[4:n]
lm(xp~0+X)

ar = arima(x, c(3,0,0), include.mean=FALSE)

# Plot ACF and PACF of fitted model to compare against true
barplot(ARMAacf(ar=c(ar$coef), lag.max=8, pacf=FALSE))
barplot(ARMAacf(ar=c(ar$coef), lag.max=8, pacf=TRUE))
