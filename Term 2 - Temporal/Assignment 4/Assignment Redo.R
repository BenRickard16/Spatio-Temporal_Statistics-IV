# Variation of GC across a gene sequence
library(astsa)
dna = bnrf1ebv # the BNRF1 gene of EBV
head(dna)

## 1=A, 2=C, 3=G, 4=T
gcTF = (dna==2)|(dna==3)
gc = gcTF*1
head(gc)

length(gc)
mean(gc)

# Linear filter to gc
gcSmooth = filter(gc, rep(1/150, 150))
tsplot(gcSmooth, col=4, lwd=2, ylim=c(0, 1))

# Fitting a HMM switching between 'typical' and 'GC-rich'
# Mean lengths of 200 and 100
# Probability of 0.6, and 0.7 of emission

# Functions required
hmmFilter = function(P, l)
  function(f, y) {
    fNew = (f %*% P) * l(y)
    fNew / sum(fNew)
  }
hmmSmoother = function(P)
  function(fp, sp) {
    fp * ((sp / (fp %*% P)) %*% t(P))
  }

# We can now use these to fit the model and plot the results.
P = matrix(c(0.995, 0.005, 0.01, 0.99), ncol=2, byrow=TRUE)
pi0 = c(0.5, 0.5)
ep = c(0.6, 0.7) # emission probabilities

## Bernoulli likelihood
l = function(y) c(dbinom(y,1,ep[1]), dbinom(y,1,ep[2]))
fpList = Reduce(hmmFilter(P, l), gc, pi0, acc=TRUE)
spList = Reduce(hmmSmoother(P), fpList, right=TRUE, acc=TRUE)
spMat = sapply(spList, cbind)
spTs = ts(spMat[2,-1])
tsplot(gcSmooth, col=4, lwd=2, ylim=c(0, 1))
lines(spTs, col=2)


# Want to forecast monthly mean carbon dioxide levels in Hawaii
tsplot(cardox, col=4, lwd=1.5, main='Monthly carbon dioxide levels')

# Fit a locally linear plus monthly seasonal effects DLM
y = cardox

library(dlm)
buildMod = function(lpar)
  dlmModPoly(2, dV=exp(lpar[1]), dW=exp(lpar[2:3])) +
    dlmModSeas(12, dV=0, dW=c(exp(lpar[4]), rep(0,10)))

opt = dlmMLE(y, parm=log(rep(0.01, 4)), build=buildMod)
opt

paste(c("dV","dW1","dW2","dWs"), "=", exp(opt$par))

mod = buildMod(opt$par)

# Using these variances, generate forecasts and intervals for the next 60 months
fit = dlmFilter(y, mod)
fore = dlmForecast(fit, 60)
pred = ts(c(tail(y, 1), fore$f),
          start=end(y), frequency=frequency(y))
upper = ts(c(tail(y, 1), fore$f + 2*sqrt(unlist(fore$Q))),
           start=end(y), frequency=frequency(y))
lower = ts(c(tail(y, 1), fore$f- 2*sqrt(unlist(fore$Q))),
           start=end(y), frequency=frequency(y))
all = ts(c(y, upper[-1]), start=start(y),
         frequency = frequency(y))
tsplot(all, xlim=c(2010, 2025), ylim=c(380, 430),
       main="Forecasts with 2SD intervals", ylab="CO2 levels")
lines(y, col=4, lwd=1.5)
lines(pred, col=2, lwd=2)
lines(upper, col=2)
lines(lower, col=2)


# Model cardiovascular mortality in LA
tsplot(cmort, col=4, lwd=1.5)

# Locally constant model with seasonal effect from 3 Fourier harmonics
# Fit model using MLE and overlay smoothed values (including and excluding seasonal effect)
# on a plot of the data
y = cmort
buildMod = function(lpar)
  dlmModPoly(1, dV=exp(lpar[1]), dW=exp(lpar[2])) +
    dlmModTrig(52, 3, dV=0, dW=exp(rep(lpar[3], 6)))

opt = dlmMLE(y, parm=rep(0,3), build=buildMod)
opt

paste(c("dV","dW","dWf"), "=", exp(opt$par))

# Smoothing and plotting
mod = buildMod(opt$par)
tsplot(y, col=4, lwd=1.5)
ss = dlmSmooth(y, mod)
so = ss$s[-1,] %*% t(mod$FF)
so = ts(so, start(y), frequency=frequency(y))
lines(ss$s[,1], col=3, lwd=2)
lines(so, col=2, lwd=2)

# What if instead we fit a locally linear model
buildMod = function(lpar)
  dlmModPoly(2, dV=exp(lpar[1]), dW=exp(lpar[2:3])) +
  dlmModTrig(52, 3, dV=0, dW=exp(rep(lpar[4], 6)))

opt2 = dlmMLE(y, parm=rep(0,4), build=buildMod)
opt2

paste(c("dV","dW1", "dW2", "dWf"), "=", exp(opt2$par))

mod = buildMod(opt2$par)
tsplot(y, col=4, lwd=1.5)
ss = dlmSmooth(y, mod)
so = ss$s[-1,] %*% t(mod$FF)
so = ts(so, start(y), frequency=frequency(y))
lines(ss$s[,1], col=3, lwd=2)
lines(so, col=2, lwd=2)
# Log-likelihood is slightly higher (expected as nested models), but visually
# the fit looks worse (not smooth) - probably not worth the additional complexity

# What if we add an AR(2) layer into model dynamics
buildMod = function(lpar)
  dlmModPoly(1, dV=exp(lpar[1]), dW=exp(lpar[2])) +
  dlmModTrig(52, 3, dV=0, dW=exp(rep(lpar[3], 6))) +
  dlmModARMA(ar=lpar[4:5], sigma2=exp(lpar[6]), dV=0)

opt3 = dlmMLE(y, parm=c(opt$par, 0.1, 0, log(0.01)),
              build=buildMod)
opt3

mod = buildMod(opt3$par)
tsplot(y, col=4, lwd=1.5)
ss = dlmSmooth(y, mod)
so = ss$s[-1,] %*% t(mod$FF)
so = ts(so, start(y), frequency=frequency(y))
lines(ss$s[,1], col=3, lwd=2)
lines(so, col=2, lwd=2)

# Minor improvement in likelihood over original model, and the fit looks fine, but the
# estimated innovation variance of the AR(2) layer is very small so not adding much