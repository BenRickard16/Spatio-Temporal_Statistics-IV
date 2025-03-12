### Spectral Methods and Digital Audio ###

## Question 1

# Spectral density of an AR(2) model
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



## Question 2

# Pure Tones
# Create a time series for middle C for 0.2 seconds
freq = 44100
times = seq(0, 0.2, by=1/freq)
midC = sin(2*pi*262*times)
plot(ts(midC, start=0, freq=freq), col=4)

# Calculating the DFT
length(midC)
spec = abs(fft(midC))
length(spec)
plot(ts(spec), col=3)

# Zooming on the first peak (since mirror image)
plot(ts(spec[1:150]), col=3)
which.max(spec[1:150])
# So single sharp peak at DFT coefficient index of 53, which is DFT coefficient 52
# SO frequency in times series time steps of
52/8821
# But time series itself has frequency of 44.1kHz, this is a frequency in Hz of
44100*52/8821

# Installing sound packahe
install.packages('sound')
library(sound)
help(package='sound')

# Turn time series representation of a sound into a sound sample object
midCsamp = as.Sample(midC)
str(midCsamp)

play(midCsamp)
saveSample(midCsamp, 'midC.wav', overwrite=TRUE)



# Instead could have generated tone
midCs <- Sine(262, 0.2)
play(midCs)

# Synthesise 2 pure notes separated by a gap
tune <- appendSample(Sine(262, 0.2), Silence(0.1, rate=freq), Sine(440, 0.3))
tune
plot(tune)

# More complicated tune
semitone <- function(h, base = 440) base*2^(h/12)

tune <- appendSample(
  Sine(semitone(10), 0.2), Silence(0.1, rate=1440),
  Sine(semitone(12), 0.2), Silence(0.1, rate=1440),
  Sine(semitone(8), 0.2), Silence(0.1, rate=1440),
  Sine(semitone(-4), 0.2), Silence(0.1, rate=1440),
  Sine(semitone(3), 0.4)
)

play(tune)
saveSample(tune, 'tune1.wav', overwrite=TRUE)


# Touche-Tone Phones
# Tone for number 1
one <- 0.5*(Sine(697, 0.2) + Sine(1209, 0.2))
plot(one)
plot(one[1:1000])

# Checking the DFT of sample
data <- sound(one)[1,]
spec <- abs(fft(data))
plot(ts(spec), col=3)
plot(ts(spec[1:500]), col=3)

# Cut-off frequency 1kHz corresponds to index
1000 * length(data) / freq

# Row index
which.max(spec[1:200])

# Column index
200 + which.max(spec[201:500])

# Correspond to actual frequencies
139 * freq / length(data)
242 * freq / length(data)









