# Set output types
region = 'US'
n = 500000
date = "2020-08-25"
options(scipen = 10)

## Load COVID tracking data
load("./manuscript version/import/covidtracking-2020-08-22.rdata")

# Assign observed data according to region
if(region == 'US') out.region <- out_daily$y_us
if(region == 'NY') out.region <- cbind(matrix(NA, nr = 5, nc = 42), out_daily$y_ny)
if(region == 'CA') out.region <- cbind(matrix(NA, nr = 5, nc = 42), out_daily$y_ca)

# Smooth out testing data
require(forecast)
out.y <- out.region
for(i in 1:dim(out.y)[1]) out.y[i,] <- as.numeric(ma(out.y[i,], order = 7))

# Source COVID pandemic evolution functions
source("./manuscript version/programs/covid_functions.r")

# Load posterior samples
load(file = paste0("./manuscript version/output/covid_pf_API_", region, "_", n, "-", date, ".rdata"))

# Calculate 95% credible intervals of marginal filtered distributions
nage <- (dim(out.samp$state)[1] - 3) / 21
require(Hmisc)
theta.wq <- array(NA, c(2, dim(out.samp$theta)[3], dim(out.samp$theta)[1]))
theta.m <- matrix(nr = dim(out.samp$theta)[1], nc = dim(out.samp$theta)[3])
out.wq <- array(NA, c(2, dim(out.samp$state)[3], dim(out.samp$state)[1]))
out.m <- matrix(nr = dim(out.samp$state)[1], nc = dim(out.samp$state)[3])
state.wq <- array(NA, c(2, dim(out.samp$state)[3], 9))
state.m <- matrix(nr = 9, nc = dim(out.samp$state)[3])
for(i in 1:dim(out.samp$state)[3]){
  # Posterior means and 95% CrI for fixed unknown parameters
  theta.wq[,i,c(1, 3:6, 12:15, 17, 19:21)] <- apply(theta2u(out.samp$theta[c(1, 3:6, 12:15, 17, 19:21),,i], 0, 1), 1, function(x) quantile(x, probs = c(0.025, .975)))
  theta.wq[,i,c(2, 7:11, 16, 18, 22:25)] <- apply(exp(out.samp$theta[c(2, 7:11, 16, 18, 22:25),,i]), 1, function(x) quantile(x, probs = c(0.025, .975)))
  theta.m[c(1, 3:6, 12:15, 17, 19:21),i] <- apply(theta2u(out.samp$theta[c(1, 3:6, 12:15, 17, 19:21),,i], 0, 1), 1, mean)
  theta.m[c(2, 7:11, 16, 18, 22:25),i] <- apply(exp(out.samp$theta[c(2, 7:11, 16, 18, 22:25),,i]), 1, mean)
  
  # Posterior means and 95% CrI for all individual states of model (by each age group)
  out.wq[,i,] <- apply(out.samp$state[,,i], 1, function(x) quantile(x, probs = c(0.025, .975)))
  out.m[,i] <- apply(out.samp$state[,,i], 1, mean)
  
  # Posterior means and 95% CrI for overall SEIR states and resource use
  state.wq[,i,1] <- quantile(apply(out.samp$state[1:nage,,i], 2, sum), probs = c(0.025, .975)) # S
  state.m[1,i] <- mean(apply(out.samp$state[1:nage,,i], 2, sum))
  state.wq[,i,2] <- quantile(apply(out.samp$state[(nage+1):(2*nage),,i], 2, sum), probs = c(0.025, .975)) # E
  state.m[2,i] <- mean(apply(out.samp$state[(nage+1):(2*nage),,i], 2, sum))
  state.wq[,i,3] <- quantile(apply(out.samp$state[c((2*nage+1):(14*nage), (15*nage+1):(18*nage)),,i], 2, sum), probs = c(0.025, .975)) # I
  state.m[3,i] <- mean(apply(out.samp$state[c((2*nage+1):(14*nage), (15*nage+1):(18*nage)),,i], 2, sum))
  state.wq[,i,4] <- quantile(apply(out.samp$state[(20*nage+1):(21*nage),,i], 2, sum), probs = c(0.025, .975)) # R
  state.m[4,i] <- mean(apply(out.samp$state[(20*nage+1):(21*nage),,i], 2, sum))
  state.wq[,i,5] <- quantile(apply(out.samp$state[(14*nage+1):(15*nage),,i], 2, sum), probs = c(0.025, .975)) # Positives
  state.m[5,i] <- mean(apply(out.samp$state[(14*nage+1):(15*nage),,i], 2, sum))
  state.wq[,i,6] <- quantile(apply(out.samp$state[(15*nage+1):(18*nage),,i], 2, sum), probs = c(0.025, .975)) # All hosp.
  state.m[6,i] <- mean(apply(out.samp$state[(15*nage+1):(18*nage),,i], 2, sum))
  state.wq[,i,7] <- quantile(apply(out.samp$state[(16*nage+1):(18*nage),,i], 2, sum), probs = c(0.025, .975)) # All ICU
  state.m[7,i] <- mean(apply(out.samp$state[(16*nage+1):(18*nage),,i], 2, sum))
  state.wq[,i,8] <- quantile(apply(out.samp$state[(17*nage+1):(18*nage),,i], 2, sum), probs = c(0.025, .975)) # Vent
  state.m[8,i] <- mean(apply(out.samp$state[(17*nage+1):(18*nage),,i], 2, sum))
  state.wq[,i,9] <- quantile(apply(out.samp$state[(19*nage+1):(20*nage),,i], 2, sum), probs = c(0.025, .975)) # New hosp.
  state.m[9,i] <- mean(apply(out.samp$state[(19*nage+1):(20*nage),,i], 2, sum))
  
  print(i)
}

# Set plotting parameters
t0 <- 1
tmax <- dim(out.samp$state)[3]
t <- seq(as.Date("1/20/2020", "%m/%d/%y") + t0, as.Date("1/20/2020", "%m/%d/%y") + tmax, "days")

# Plot dynamic states over time
png(filename = paste0("./manuscript version/graphs/covid_pf_", region, "_", n,"_dynamicStates-",Sys.Date(),".png"),
    width = 800, height = 583)

# Plot social-distancing parameter over time
par(mfrow=c(3,1))
plot(t, out.m[21*nage + 1,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(xi),
     ylim = c(min(out.m[21*nage + 1,t0:tmax], out.wq[1,t0:tmax,21*nage + 1]),
              max(out.m[21*nage + 1,t0:tmax], out.wq[2,t0:tmax,21*nage + 1])))
lines(t, out.wq[1,t0:tmax,21*nage + 1], lty = 'dashed')
lines(t, out.wq[2,t0:tmax,21*nage + 1], lty = 'dashed')

# Plot testing rates over time
plot(t, out.m[21*nage + 2,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(w[A]),
     ylim = c(min(out.m[21*nage + 2,t0:tmax], out.wq[1,t0:tmax,21*nage + 2]),
              max(out.m[21*nage + 2,t0:tmax], out.wq[2,t0:tmax,21*nage + 2])))
lines(t, out.wq[1,t0:tmax,21*nage + 2], lty = 'dashed')
lines(t, out.wq[2,t0:tmax,21*nage + 2], lty = 'dashed')

plot(t, out.m[21*nage + 3,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(w[I]),
     ylim = c(min(out.m[21*nage + 3,t0:tmax], out.wq[1,t0:tmax,21*nage + 3]),
              max(out.m[21*nage + 3,t0:tmax], out.wq[2,t0:tmax,21*nage + 3])))
lines(t, out.wq[1,t0:tmax,21*nage + 3], lty = 'dashed')
lines(t, out.wq[2,t0:tmax,21*nage + 3], lty = 'dashed')

dev.off()

# Plot observed versus true versus estimated resource use and deaths
dat.sird.mean <- rbind(data.frame(Time = t, N = state.m[5,t0:tmax], state = "Positive tests"),
                       data.frame(Time = t, N = state.m[6,t0:tmax], state = 'In hospital'),
                       data.frame(Time = t, N = state.m[7,t0:tmax], state = 'In ICU'),
                       data.frame(Time = t, N = state.m[8,t0:tmax], state = 'On ventilator'),
                       data.frame(Time = t, N = state.m[9,t0:tmax], state = 'Died'))

dat.sird.lb <- rbind(data.frame(Time = t, N = state.wq[1,t0:tmax,5], state = 'Positive tests'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,6], state = 'In hospital'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,7], state = 'In ICU'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,8], state = 'On ventilator'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,9], state = 'Died'))

dat.sird.ub <- rbind(data.frame(Time = t, N = state.wq[2,t0:tmax,5], state = 'Positive tests'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,6], state = 'In hospital'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,7], state = 'In ICU'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,8], state = 'On ventilator'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,9], state = 'Died'))

dat.sird.observed <- rbind(data.frame(Time = t[-1], N = out.region[1,t0:(tmax-1)], state = 'Positive tests'),
                           data.frame(Time = t[-1], N = out.region[2,t0:(tmax-1)], state = 'In hospital'),
                           data.frame(Time = t[-1], N = out.region[3,t0:(tmax-1)], state = 'In ICU'),
                           data.frame(Time = t[-1], N = out.region[4,t0:(tmax-1)], state = 'On ventilator'),
                           data.frame(Time = t[-1], N = out.region[5,t0:(tmax-1)], state = 'Died'))

dat.sird <- rbind(data.frame(dat.sird.mean, Type = "Posterior Mean", Measure = "Posterior Mean"),
                  data.frame(dat.sird.lb, Type = "95% CI Lower Bound", Measure = "95% CI"),
                  data.frame(dat.sird.ub, Type = "95% CI Upper Bound", Measure = "95% CI"),
                  data.frame(dat.sird.observed, Type = "Observed", Measure = "Observed"))
dat.sird$Type <- factor(dat.sird$Type, levels = c("Observed", "Posterior Mean", "95% CI Lower Bound", "95% CI Upper Bound"))
dat.sird$Measure <- factor(dat.sird$Measure, levels = c("Observed", "Posterior Mean", "95% CI"))

require(ggplot2)
p1 <- ggplot(data = dat.sird, aes(x=Time,y=N)) + geom_line(aes(group=Type, col=Measure))
g <- p1 + facet_wrap(~state, nrow = 2, scales = 'free_y') + 
  labs(x = "Time", y = "Number of Persons") + scale_x_date(date_labels = "%b",date_breaks = "1 month") +
  theme_bw() + theme(strip.text.x=element_text(size=16),strip.text.y=element_text(size=16),
                     plot.title=element_text(size=20),
                     axis.title=element_text(size = 15), axis.text=element_text(size=12),
                     legend.title=element_text(size=15),legend.text=element_text(size=12))
ggsave(file = paste0("./manuscript version/graphs/covid_pf_", region, "_", n,"_ResourceObserved-",Sys.Date(),".png"), g, width = 12, height = 7)

# Plot SEIR curves
dat.sird.mean <- rbind(data.frame(Time = t, N = state.m[1,t0:tmax], state = 'Susceptible'),
                       data.frame(Time = t, N = state.m[2,t0:tmax], state = 'Exposed'),
                       data.frame(Time = t, N = state.m[3,t0:tmax], state = 'Infected'),
                       data.frame(Time = t, N = state.m[4,t0:tmax], state = 'Recovered'))

dat.sird.lb <- rbind(data.frame(Time = t, N = state.wq[1,t0:tmax,1], state = 'Susceptible'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,2], state = 'Exposed'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,3], state = 'Infected'),
                     data.frame(Time = t, N = state.wq[1,t0:tmax,4], state = 'Recovered'))

dat.sird.ub <- rbind(data.frame(Time = t, N = state.wq[2,t0:tmax,1], state = 'Susceptible'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,2], state = 'Exposed'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,3], state = 'Infected'),
                     data.frame(Time = t, N = state.wq[2,t0:tmax,4], state = 'Recovered'))

dat.sird <- rbind(data.frame(dat.sird.mean, Type = "Posterior Mean", Measure = "Posterior Mean"),
                  data.frame(dat.sird.lb, Type = "95% CI Lower Bbound", Measure = "95% CI"),
                  data.frame(dat.sird.ub, Type = "95% CI Upper Bound", Measure = "95% CI"))

require(ggplot2)
p1 <- ggplot(data = dat.sird, aes(x=Time,y=N)) + geom_line(aes(group=Type, col=Measure))
g <- p1 + facet_wrap(~state, nrow = 2, scales = 'free_y') + 
  labs(x = "Time", y = "Number of Persons") + scale_x_date(date_labels = "%b",date_breaks = "1 month") +
  theme_bw() + theme(strip.text.x=element_text(size=16),strip.text.y=element_text(size=16),
                     plot.title=element_text(size=20),
                     axis.title=element_text(size = 15), axis.text=element_text(size=12),
                     legend.title=element_text(size=15),legend.text=element_text(size=12))
ggsave(file = paste0("./manuscript version/graphs/covid_pf_", region, "_", n,"_SEIR-",Sys.Date(),".png"), g, width = 10, height = 6)

# Plot 95% CIs of marginal filtered distributions of unknown parameters
pdf(file = paste0("./manuscript version/graphs/covid_pf_", region, "_", n, "_params-",Sys.Date(),".pdf"), width = 8, height = 10)

# Testing parameters
par(mfrow=c(2,1))
plot(t, theta.m[1,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(p[alpha]),
     ylim = c(min(theta.m[1,t0:tmax], theta.wq[1,t0:tmax,1]),
              max(theta.m[1,t0:tmax], theta.wq[2,t0:tmax,1])))
lines(t, theta.wq[1,t0:tmax,1], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,1], lty = 'dashed')

plot(t, theta.m[2,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(sigma[t]),
     ylim = c(min(theta.m[2,t0:tmax], theta.wq[1,t0:tmax,2]),
              max(theta.m[2,t0:tmax], theta.wq[2,t0:tmax,2])))
lines(t, theta.wq[1,t0:tmax,2], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,2], lty = 'dashed')

# Proportion hospitalized
par(mfrow = c(4,1))
plot(t, theta.m[3,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(rho[1]),
     ylim = c(min(theta.m[3,t0:tmax], theta.wq[1,t0:tmax,3]),
              max(theta.m[3,t0:tmax], theta.wq[2,t0:tmax,3])))
lines(t, theta.wq[1,t0:tmax,3], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,3], lty = 'dashed')

plot(t, theta.m[4,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(rho[2]),
     ylim = c(min(theta.m[4,t0:tmax], theta.wq[1,t0:tmax,4]),
              max(theta.m[4,t0:tmax], theta.wq[2,t0:tmax,4])))
lines(t, theta.wq[1,t0:tmax,4], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,4], lty = 'dashed')

plot(t, theta.m[5,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(rho[3]),
     ylim = c(min(theta.m[5,t0:tmax], theta.wq[1,t0:tmax,5]),
              max(theta.m[5,t0:tmax], theta.wq[2,t0:tmax,5])))
lines(t, theta.wq[1,t0:tmax,5], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,5], lty = 'dashed')

plot(t, theta.m[6,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(rho[4]),
     ylim = c(min(theta.m[6,t0:tmax], theta.wq[1,t0:tmax,6]),
              max(theta.m[6,t0:tmax], theta.wq[2,t0:tmax,6])))
lines(t, theta.wq[1,t0:tmax,6], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,6], lty = 'dashed')

# Hospitalized variance and rate
par(mfrow = c(2,1))
plot(t, theta.m[7,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(sigma[h]),
     ylim = c(min(theta.m[7,t0:tmax], theta.wq[1,t0:tmax,7]),
              max(theta.m[7,t0:tmax], theta.wq[2,t0:tmax,7])))
lines(t, theta.wq[1,t0:tmax,7], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,7], lty = 'dashed')

plot(t, theta.m[8,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(delta),
     ylim = c(min(theta.m[8,t0:tmax], theta.wq[1,t0:tmax,8]),
              max(theta.m[8,t0:tmax], theta.wq[2,t0:tmax,8])))
lines(t, theta.wq[1,t0:tmax,8], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,8], lty = 'dashed')

# Recovery rates from hospital
par(mfrow = c(3,1))
plot(t, theta.m[9,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(psi[h]),
     ylim = c(min(theta.m[9,t0:tmax], theta.wq[1,t0:tmax,9]),
              max(theta.m[9,t0:tmax], theta.wq[2,t0:tmax,9])))
lines(t, theta.wq[1,t0:tmax,9], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,9], lty = 'dashed')

plot(t, theta.m[10,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(psi[c]),
     ylim = c(min(theta.m[10,t0:tmax], theta.wq[1,t0:tmax,10]),
              max(theta.m[10,t0:tmax], theta.wq[2,t0:tmax,10])))
lines(t, theta.wq[1,t0:tmax,10], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,10], lty = 'dashed')

plot(t, theta.m[11,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(psi[v]),
     ylim = c(min(theta.m[11,t0:tmax], theta.wq[1,t0:tmax,11]),
              max(theta.m[11,t0:tmax], theta.wq[2,t0:tmax,11])))
lines(t, theta.wq[1,t0:tmax,11], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,11], lty = 'dashed')

# Proportion in ICU among hospitalized
par(mfrow = c(4,1))
plot(t, theta.m[12,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(c[1]),
     ylim = c(min(theta.m[12,t0:tmax], theta.wq[1,t0:tmax,12]),
              max(theta.m[12,t0:tmax], theta.wq[2,t0:tmax,12])))
lines(t, theta.wq[1,t0:tmax,12], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,12], lty = 'dashed')

plot(t, theta.m[13,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(c[2]),
     ylim = c(min(theta.m[13,t0:tmax], theta.wq[1,t0:tmax,13]),
              max(theta.m[13,t0:tmax], theta.wq[2,t0:tmax,13])))
lines(t, theta.wq[1,t0:tmax,13], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,13], lty = 'dashed')

plot(t, theta.m[14,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(c[3]),
     ylim = c(min(theta.m[14,t0:tmax], theta.wq[1,t0:tmax,14]),
              max(theta.m[14,t0:tmax], theta.wq[2,t0:tmax,14])))
lines(t, theta.wq[1,t0:tmax,14], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,14], lty = 'dashed')

plot(t, theta.m[15,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(c[4]),
     ylim = c(min(theta.m[15,t0:tmax], theta.wq[1,t0:tmax,15]),
              max(theta.m[15,t0:tmax], theta.wq[2,t0:tmax,15])))
lines(t, theta.wq[1,t0:tmax,15], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,15], lty = 'dashed')

# ICU variance
par(mfrow=c(1,1))
plot(t, theta.m[16,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(sigma[i]),
     ylim = c(min(theta.m[16,t0:tmax], theta.wq[1,t0:tmax,16]),
              max(theta.m[16,t0:tmax], theta.wq[2,t0:tmax,16])))
lines(t, theta.wq[1,t0:tmax,16], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,16], lty = 'dashed')

# Proportion ventilated and ventilation variance
par(mfrow=c(2,1))
plot(t, theta.m[17,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(v),
     ylim = c(min(theta.m[17,t0:tmax], theta.wq[1,t0:tmax,17]),
              max(theta.m[17,t0:tmax], theta.wq[2,t0:tmax,17])))
lines(t, theta.wq[1,t0:tmax,17], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,17], lty = 'dashed')

plot(t, theta.m[18,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(sigma[v]),
     ylim = c(min(theta.m[18,t0:tmax], theta.wq[1,t0:tmax,18]),
              max(theta.m[18,t0:tmax], theta.wq[2,t0:tmax,18])))
lines(t, theta.wq[1,t0:tmax,18], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,18], lty = 'dashed')

# Mortality weights
par(mfrow = c(3,1))
plot(t, theta.m[19,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(m[h]),
     ylim = c(min(theta.m[19,t0:tmax], theta.wq[1,t0:tmax,19]),
              max(theta.m[19,t0:tmax], theta.wq[2,t0:tmax,19])))
lines(t, theta.wq[1,t0:tmax,19], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,19], lty = 'dashed')

plot(t, theta.m[20,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(m[c]),
     ylim = c(min(theta.m[20,t0:tmax], theta.wq[1,t0:tmax,20]),
              max(theta.m[20,t0:tmax], theta.wq[2,t0:tmax,20])))
lines(t, theta.wq[1,t0:tmax,20], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,20], lty = 'dashed')

plot(t, theta.m[21,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(m[v]),
     ylim = c(min(theta.m[21,t0:tmax], theta.wq[1,t0:tmax,21]),
              max(theta.m[21,t0:tmax], theta.wq[2,t0:tmax,21])))
lines(t, theta.wq[1,t0:tmax,21], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,21], lty = 'dashed')

# Mortality rates
par(mfrow = c(3,1))
plot(t, theta.m[22,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(mu[h]),
     ylim = c(min(theta.m[22,t0:tmax], theta.wq[1,t0:tmax,22]),
              max(theta.m[22,t0:tmax], theta.wq[2,t0:tmax,22])))
lines(t, theta.wq[1,t0:tmax,22], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,22], lty = 'dashed')

plot(t, theta.m[23,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(mu[c]),
     ylim = c(min(theta.m[23,t0:tmax], theta.wq[1,t0:tmax,23]),
              max(theta.m[23,t0:tmax], theta.wq[2,t0:tmax,23])))
lines(t, theta.wq[1,t0:tmax,23], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,23], lty = 'dashed')

plot(t, theta.m[24,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(mu[v]),
     ylim = c(min(theta.m[24,t0:tmax], theta.wq[1,t0:tmax,24]),
              max(theta.m[24,t0:tmax], theta.wq[2,t0:tmax,24])))
lines(t, theta.wq[1,t0:tmax,24], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,24], lty = 'dashed')

# Death variance
par(mfrow=c(1,1))
plot(t, theta.m[25,t0:tmax], type = 'l', xlab = expression(t), ylab = expression(sigma[d]),
     ylim = c(min(theta.m[25,t0:tmax], theta.wq[1,t0:tmax,25]),
              max(theta.m[25,t0:tmax], theta.wq[2,t0:tmax,25])))
lines(t, theta.wq[1,t0:tmax,25], lty = 'dashed')
lines(t, theta.wq[2,t0:tmax,25], lty = 'dashed')

dev.off()

# Examine trajectories of dynamic particle states and filtered resource utilization

# Load sampled particle traces
load(file = paste0("./manuscript version/output/covid_pf_API_trace_", region, "_", n, "-", date, ".rdata"))

# How many particles were sampled
n.samp <- dim(out.trace)[2]

# Set plot variables
pos.mean <- 5:9
trace <- list((14*nage + 1):(15*nage), (15*nage + 1):(18*nage), (16*nage + 1):(18*nage),
              (17*nage + 1):(18*nage), (19*nage + 1):(20*nage))
ylab = c("Daily positive tests", "In hospital", "In ICU", "On ventilator", "Died")

# Find particle trajectory that minimizes the sum of squared differences between predicted and observed data
ssmin <- Inf
for(i in 1:n.samp){
  ss <- 0
  for(j in 1:dim(out.y)[1]) ss <- ss + sum((out.y[j,] - apply(out.trace[trace[[j]],i,-1], 2, sum))^2, na.rm=TRUE)
  if(ss < ssmin){
    part.min <- i
    ssmin <- ss
  }
}

# Plot trajectories -- resource use
png(filename = paste0("./manuscript version/graphs/covid_pf_trace_", region, "_", n, "-ResourceObserved-", Sys.Date(), ".png"),
    width = 1000, height = 583)

par(mfrow=c(2,3), mar = c(5,5,4,2)+0.1)

for(j in 1:5){
  trace.j <- apply(out.trace[trace[[j]],,], c(2,3), sum)
  plot(t, state.m[pos.mean[j],t0:tmax], ylim = c(min(state.m[pos.mean[j],], out.region[j,], trace.j, na.rm=TRUE),
                                                 max(state.m[pos.mean[j],], out.region[j,], trace.j, na.rm = TRUE)),
       type = 'l', col='white', xlab = expression(t), ylab = ylab[j], cex.lab = 2)
  if(j == 2) title(main = "Particle Trace of Resource Use", cex.main = 1.8)
  for(i in 1:n.samp) lines(t, trace.j[i,t0:tmax], col = 'gray', lwd = .5)
  lines(t, state.m[pos.mean[j],t0:tmax], lty = 2)
  lines(t[-1], out.region[j,t0:(tmax-1)], col = 2)
  lines(t, apply(apply(out.trace[trace[[j]],,t0:tmax], c(2, 3), sum), 2, mean))
  lines(t, apply(out.trace[trace[[j]],part.min,t0:tmax], 2, sum), lwd = 2)
  if(j == 1) legend("topright", legend = c("Trace", "Filtered mean", "Observed", "Mean trace", "Min trace"),
                    lty = c(1,2,1,1,1), col = c('gray', 'black', 'red', 'black', 'black'), lwd = c(.5, 1, 1, 1, 2), cex = 1.25)
  
}

dev.off()

# Plot trajectories -- dynamic states
png(filename = paste0("./manuscript version/graphs/covid_pf_trace_", region, "_", n, "-dynamicStates-", Sys.Date(), ".png"),
    width = 800, height = 700)

par(mfrow=c(3,1), mar = c(5,5,4,2)+0.1)

# Set plot variables
param <- c(21*nage + 1, 21*nage + 2, 21*nage + 3)
ylab = quote(c(expression(xi), expression(w[A]), expression(w[I]))[j])

for(j in 1:3){
  trace.j <- out.trace[param[j],,]
  plot(t, out.m[param[j],t0:tmax], ylim = c(min(out.m[param[j],], trace.j, na.rm=TRUE),
                                            max(out.m[param[j],], trace.j, na.rm = TRUE)),
       type = 'l', col='white', xlab = expression(t), ylab = eval(ylab), cex.lab = 2)
  if(j == 1) title(main = "Particle Trace of Dynamic States", cex.main = 1.8)
  for(i in 1:n.samp) lines(t, trace.j[i,t0:tmax], col = 'gray', lwd = .5)
  lines(t, out.m[param[j],t0:tmax], lty = 2)
  lines(t,  apply(out.trace[param[j],,t0:tmax], 2, mean))
  lines(t, out.trace[param[j],part.min,t0:tmax], lwd = 2)
  if(j == 1) legend("topright", legend = c("Trace", "Filtered mean", "Mean trace", "Min trace"),
                    lty = c(1,2,1,1), col = c('gray', 'black', 'black', 'black'), lwd = c(.5, 1, 1, 2), cex = 1.15)
}

dev.off()

# Run model using mean trace and posterior means for fixed parameters

# Set values of known fixed parameters
# Population size
if(region == 'US') P <- 330000000 * c(.25, .4, .19, .16)
if(region == 'NY') P <- 19542209 * c(.23, .4, .20, .17)
if(region == 'CA') P <- 39557045 * c(.25, .42, .19, .14)

# Contact matrices
M <- cbind(c(9.76, 3.77, 1.51, 0.6), c(3.77, 9.43, 3.05, 0.7), c(1.51, 3.05, 2.96, 0.76), c(0.6, 0.7, .76, 1.25))
Mtilde <- cbind(c(2.04, 1.56, 0.5, 0.38), c(1.56, 1.51, 0.45, 0.24), c(0.5, 0.45, 1.04, 0.19), c(0.38, 0.24, 0.19, 0.64))

# Proportion of asymptomatic cases
zeta <- rep(.179, 4)

# Proportion of mild cases out of those with symptoms
theta.known <- c(0.973, 0.90, 0.648, 0.466)

# Transmission parameters
beta <- 0.0493 # transmission rate per contact with infectious individual
kappa_M <- 0.5 # relative infectivity of a mild case relative to severe
kappa_A <- 0.5 # relative infectivity of an asymptomatic case relative to severe
sigma <- 1 / 5.2 # inverse of the average incubation period
gamma <- 1 / 4.6 # inverse of the average recovery period

# Behavioral parameters
q <- 0.05 # proportion of severe cases that self-isolate immediately upon symptom onset
f_A <- 0.05 # baseline self-isolation factor for mild cases
tau_A <- 1 / 2 # inverse of the average time until self-isolation for mild cases
f_I <- 0.8 # baseline self-isolation factor for severe cases
tau_I <- 1 / 1 # inverse of the average time until self-isolation for severe cases
f_T <- 0.7 # proportion of cases that self-isolate while awaiting test results
phi_x <- 0.95 # autocorrelation parameter for local transmission / social-distancing adjustment
sigma_x <- 0.25 # standard deviation parameter for local transmission / social-distancing adjustment

# Testing parameters
w_s <- 0.1 # proportion of asymptomatic cases that get tested
r_A <- 1 / 3 # inverse of the average time from symptom onset to test for mild cases
r_I <- 1 / 2 # inverse of the average time from symptom onset to test for severe cases
phi <- 1 / 2.25885149849582 # inverse of the average turnaround time for testing
omega <- 0.972360338537735 # sensitivity of testing
phi_A <- 0.995 # autocorrelation parameter for proportion of mild cases tested
phi_I <- 0.99 # autocorrelation parameter for proportion of severe cases tested
sigma_w <- 0.01 # standard deviation parameter for mild / severe cases tested

# Set unknown fixed parameters at posterior means
if(region == 'US') alpha <- theta.m[1, dim(theta.m)[2]]*1122897.86964862 # daily throughput cap
if(region == 'NY')  alpha <- theta.m[1, dim(theta.m)[2]]*66550.0668668088 # daily throughput cap
if(region == 'CA') alpha <- theta.m[1, dim(theta.m)[2]]*135170.166670578 # daily throughput cap
rho <- theta.m[3:6, dim(theta.m)[2]]
delta <- theta.m[8, dim(theta.m)[2]]
psi_h <- theta.m[9, dim(theta.m)[2]]
psi_c <- theta.m[10, dim(theta.m)[2]]
psi_v <- theta.m[11, dim(theta.m)[2]]
c <- theta.m[12:15, dim(theta.m)[2]]
v <- theta.m[17, dim(theta.m)[2]]
m_h <- theta.m[19, dim(theta.m)[2]]
m_c <- theta.m[20, dim(theta.m)[2]]
m_v <- theta.m[21, dim(theta.m)[2]]
mu_h <- theta.m[22, dim(theta.m)[2]]
mu_c <- theta.m[23, dim(theta.m)[2]]
mu_v <- theta.m[24, dim(theta.m)[2]]

# Set dynamic states to smoothed means
# xi.fixed <- apply(out.trace[param[1],,], 2, mean)
# w_A.fixed <- apply(out.trace[param[2],,], 2, mean)
# w_I.fixed <- apply(out.trace[param[3],,], 2, mean)
xi.fixed <- out.trace[param[1],part.min,]
w_A.fixed <- out.trace[param[2],part.min,]
w_I.fixed <- out.trace[param[3],part.min,]

# Flexible functions for particle filters
revo = function(x, theta, random = TRUE, loglik = c('t', 'h', 'i', 'v', 'd'),
                xi.fixed = NULL, w_A.fixed = NULL, w_I.fixed = NULL){
  nage <- (length(x) - 3) / 21
  x.list <- covid_prop(P = P, M = M, Mtilde = Mtilde,
                       S = x[1:nage],
                       E = x[(nage + 1):(2*nage)],
                       Isev = x[(2*nage + 1):(3*nage)],
                       Isevhosp = x[(3*nage + 1):(4*nage)],
                       Iseviso = x[(4*nage + 1):(5*nage)],
                       Isevhospiso = x[(5*nage + 1):(6*nage)],
                       Aasym = x[(6*nage + 1):(7*nage)],
                       Aasymiso = x[(7*nage + 1):(8*nage)],
                       Amild = x[(8*nage + 1):(9*nage)],
                       Amildiso = x[(9*nage + 1):(10*nage)],
                       Tasym = x[(10*nage + 1):(11*nage)],
                       Tmild = x[(11*nage + 1):(12*nage)],
                       Tsev = x[(12*nage + 1):(13*nage)],
                       Tsevhosp = x[(13*nage + 1):(14*nage)],
                       Hfloor = x[(15*nage + 1):(16*nage)],
                       Hicu = x[(16*nage + 1):(17*nage)],
                       Hvent = x[(17*nage + 1):(18*nage)],
                       D = x[(19*nage + 1):(20*nage)],
                       R = x[(20*nage + 1):(21*nage)],
                       xi = x[21*nage + 1],
                       w_A = x[21*nage + 2],
                       w_I = x[21*nage + 3],
                       
                       zeta = zeta,
                       theta = theta.known,
                       rho = if('h' %in% loglik) theta2u(theta[3:(nage+2)], 0, 1) else rho,
                       c = if('i' %in% loglik) theta2u(theta[(nage + 8):(nage + 11)], 0, 1) else c,
                       beta = beta,
                       kappa_A = kappa_A,
                       kappa_M = kappa_M,
                       sigma = sigma,
                       gamma = gamma,
                       q = q,
                       f_A = f_A,
                       tau_A = tau_A,            
                       f_I = f_I,
                       tau_I = tau_I,
                       f_T = f_T,
                       phi_x = phi_x,
                       sigma_x = sigma_x,
                       phi_A = phi_A,
                       phi_I = phi_I,
                       sigma_w = sigma_w,
                       w_s = w_s,
                       r_A = r_A,
                       r_I = r_I,
                       phi = phi,
                       omega = omega,
                       alpha = if('t' %in% loglik) theta2u(theta[1], 0, 1)*alpha else alpha,
                       delta = if('h' %in% loglik) exp(theta[nage + 4]) else delta,
                       v = if('v' %in% loglik) theta2u(theta[2*nage + 9], 0, 1) else v,
                       m_h = if('d' %in% loglik) theta2u(theta[2*nage + 11], 0, 1) else m_h,
                       m_c = if('d' %in% loglik) theta2u(theta[2*nage + 12], 0, 1) else m_c,
                       m_v = if('d' %in% loglik) theta2u(theta[2*nage + 13], 0, 1) else m_v,
                       mu_h = if('d' %in% loglik) exp(theta[2*nage + 14]) else mu_h,
                       mu_c = if('d' %in% loglik) exp(theta[2*nage + 15]) else mu_c,
                       mu_v = if('d' %in% loglik) exp(theta[2*nage + 16]) else mu_v,
                       psi_h = if('h' %in% loglik) exp(theta[nage + 5]) else psi_h,
                       psi_c = if('h' %in% loglik) exp(theta[nage + 6]) else psi_c,
                       psi_v = if('h' %in% loglik) exp(theta[nage + 7]) else psi_v,
                       
                       random = random,
                       xi.fixed = xi.fixed,
                       w_A.fixed = w_A.fixed,
                       w_I.fixed = w_I.fixed
  )
  return(c(x.list$S, x.list$E, x.list$Isev, x.list$Isevhosp, x.list$Iseviso, x.list$Isevhospiso,
           x.list$Aasym, x.list$Aasymiso, x.list$Amild, x.list$Amildiso,
           x.list$Tasym, x.list$Tmild, x.list$Tsev, x.list$Tsevhosp, x.list$Tpos,
           x.list$Hfloor, x.list$Hicu, x.list$Hvent, x.list$Hnew,
           x.list$D, x.list$R,
           x.list$xi, x.list$w_A, x.list$w_I))
}

dllik = function(y, x, theta, loglik = c('t', 'h', 'i', 'v','d')){
  nage <- (length(x) - 3) / 21
  return(covid_dllik(y = y, loglik = loglik,
                     Tpos = x[(14*nage + 1):(15*nage)],
                     Hfloor = x[(15*nage + 1):(16*nage)],
                     Hicu = x[(16*nage + 1):(17*nage)],
                     Hvent = x[(17*nage + 1):(18*nage)],
                     D = x[(19*nage + 1):(20*nage)],
                     sigma_pos = if('t' %in% loglik) exp(theta[2]) else sigma_pos,
                     sigma_hosp = if('h' %in% loglik) exp(theta[nage+3]) else sigma_hosp,
                     sigma_icu = if('i' %in% loglik) exp(theta[2*nage+8]) else sigma_icu,
                     sigma_vent = if('v' %in% loglik) exp(theta[2*nage+10]) else sigma_vent,
                     sigma_d = if('d' %in% loglik) exp(theta[2*nage+17]) else sigma_d
  ))
}

# Run APF for known states and known fixed parameters (trace model)
n.trace = 1
rprior_a <- function() rprior(loglik = c(), P=P, xi.fixed = xi.fixed, w_A.fixed = w_A.fixed, w_I.fixed = w_I.fixed)
revo_a <- function(x) revo(x, theta = NULL, loglik = c(), xi.fixed = xi.fixed, w_A.fixed = w_A.fixed, w_I.fixed = w_I.fixed)
dllik_a <- function(y, x) dllik(y, x, theta = NULL, loglik = c())
pstate <- function(x) revo(x, theta = NULL, random = FALSE, loglik = c(), xi.fixed = xi.fixed, w_A.fixed = w_A.fixed, w_I.fixed = w_I.fixed)
source("./manuscript version/programs/apf.r")
out.model.trace = apf(out.y, dllik_a, pstate, revo_a, rprior_a, n.trace)

# Plot observed versus true versus estimated resource use and deaths
dat.sird.mean <- rbind(data.frame(Time = t, N = apply(out.model.trace$state[(14*nage + 1):(15*nage),1,t0:tmax], 2, sum),
                                  state = "Positive tests"),
                       data.frame(Time = t, N = apply(out.model.trace$state[(15*nage + 1):(18*nage),1,t0:tmax], 2, sum),
                                  state = 'In hospital'),
                       data.frame(Time = t, N = apply(out.model.trace$state[(16*nage + 1):(18*nage),1,t0:tmax], 2, sum),
                                  state = 'In ICU'),
                       data.frame(Time = t, N = apply(out.model.trace$state[(17*nage + 1):(18*nage),1,t0:tmax], 2, sum),
                                  state = 'On ventilator'),
                       data.frame(Time = t, N = apply(out.model.trace$state[(19*nage + 1):(20*nage),1,t0:tmax], 2, sum),
                                  state = 'Died'))

dat.sird.observed <- rbind(data.frame(Time = t[-1], N = out.region[1,t0:(tmax-1)], state = 'Positive tests'),
                           data.frame(Time = t[-1], N = out.region[2,t0:(tmax-1)], state = 'In hospital'),
                           data.frame(Time = t[-1], N = out.region[3,t0:(tmax-1)], state = 'In ICU'),
                           data.frame(Time = t[-1], N = out.region[4,t0:(tmax-1)], state = 'On ventilator'),
                           data.frame(Time = t[-1], N = out.region[5,t0:(tmax-1)], state = 'Died'))

dat.sird <- rbind(data.frame(dat.sird.mean, Type = "Predicted", Measure = "Predicted"),
                  data.frame(dat.sird.observed, Type = "Observed", Measure = "Observed"))
dat.sird$Type <- factor(dat.sird$Type, levels = c("Observed", "Predicted"))
dat.sird$Measure <- factor(dat.sird$Measure, levels = c("Observed","Predicted"))

require(ggplot2)
p1 <- ggplot(data = dat.sird, aes(x=Time,y=N)) + geom_line(aes(group=Type, col=Measure))
g <- p1 + facet_wrap(~state, nrow = 2, scales = 'free_y') + 
  labs(x = "Time", y = "Number of Persons") + scale_x_date(date_labels = "%b",date_breaks = "1 month") +
  theme_bw() + theme(strip.text.x=element_text(size=16),strip.text.y=element_text(size=16),
                     plot.title=element_text(size=20),
                     axis.title=element_text(size = 15), axis.text=element_text(size=12),
                     legend.title=element_text(size=15),legend.text=element_text(size=12))
ggsave(file = paste0("./manuscript version/graphs/covid_pf_", region, "_", n, "_observedTrace-",Sys.Date(),".png"), g, width = 12, height = 7)

# Output smoothed dynamic states and posterior means of unknown fixed parameters
states.out <- data.frame(xi = xi.fixed, w_A = w_A.fixed, w_I = w_I.fixed)
theta.out <- theta.m[c(1, 3:6, 8:15, 17, 19:24),dim(theta.m)[2]]
names(theta.out) <- c("p_alpha", paste0("rho", 1:4), "delta", paste0("psi", c("_h", "_c", "_v")), paste0("c", 1:4), "v",
                      paste0("m", c("_h", "_c", "_v")), paste0("mu", c("_h", "_c", "_v")))
write.csv(states.out, file = paste0("./manuscript version/output/covid_pf_", region, "_", n, "_states-",Sys.Date(),".csv"),
          row.names = FALSE)
write.csv(data.frame(theta.out), file = paste0("./manuscript version/output/covid_pf_", region, "_", n, "_params-",Sys.Date(),".csv"))