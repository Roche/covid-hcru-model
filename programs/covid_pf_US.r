## Load COVID tracking data
load("./import/covidtracking-2020-08-22.rdata")

## Source COVID pandemic evolution functions
source("./programs/covid_functions.r")

# Set output labels
region = 'US'
n = 500000
options(scipen = 10)

# Assign observed data according to region
if(region == 'US') out.region <- out_daily$y_us
if(region == 'NY') out.region <- cbind(matrix(NA, nr = 5, nc = 42), out_daily$y_ny)
if(region == 'CA') out.region <- cbind(matrix(NA, nr = 5, nc = 42), out_daily$y_ca)

# Smooth out data (7-day moving average)
require(forecast)
out.y <- out.region
for(i in 1:dim(out.y)[1]) out.y[i,] <- as.numeric(ma(out.y[i,], order = 7))

## Set values of 'known' fixed parameters

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
if(region == 'US') alpha <- 1122897.86964862 # daily throughput cap
if(region == 'NY')  alpha <- 66550.0668668088 # daily throughput cap
if(region == 'CA') alpha <- 135170.166670578 # daily throughput cap
phi_A <- 0.995 # autocorrelation parameter for proportion of mild cases tested
phi_I <- 0.99 # autocorrelation parameter for proportion of severe cases tested
sigma_w <- 0.1 # standard deviation parameter for mild / severe cases tested

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

## Run KDPF for dynamic transmission states and unknown hospitalization parameters
rprior_d <- function() rprior(loglik = c('t', 'h', 'i', 'v', 'd'), P=P)
revo_d <- function(x, theta) revo(x, theta, loglik = c('t', 'h', 'i', 'v', 'd'))
dllik_d <- function(y, x, theta) dllik(y, x, theta, loglik = c('t', 'h', 'i', 'v', 'd'))
pstate <- function(x, theta) revo(x, theta, random = FALSE, loglik = c('t', 'h', 'i', 'v', 'd'))
source("./programs/kd_pf.r")
set.seed(56)
out = kd_pf(out.y + 1, dllik_d, pstate, revo_d, rprior_d, n, print = TRUE, 
            method = "stratified", nonuniformity = "ess", threshold = 0.8)

# Sample from posterior distribution of states and unknown parameters
n.samp <- 100
out.samp <- list(state = array(NA, dim = c(dim(out$state)[1], n.samp, dim(out$state)[3])),
                 theta = array(NA, dim = c(dim(out$theta)[1], n.samp, dim(out$theta)[3])))
for(i in 1:dim(out$state)[3]){
  ii <- sample(1:dim(out$state)[2], n.samp, out$weight[,i], replace = TRUE)
  out.samp$state[,,i] <- out$state[,ii,i]
  out.samp$theta[,,i] <- out$theta[,ii,i]
  print(i)
}
save(out.samp, file = paste0("./output/covid_pf_API_", region, "_", n, "-", Sys.Date(), ".rdata"))

# Find the particle paths from filtered distribution
nage <- (dim(out$state)[1] - 3) / 21
n.days <- dim(out$state)[3] - 1
n.keep <- dim(out$state)[2]
out.trace <- array(NA, dim = c(dim(out$state)[1], n.samp, dim(out$state)[3]))
ii <- sample(1:n.keep, n.samp, out$weight[,dim(out$state)[3]], replace = TRUE)
for(part.start in 1:n.samp)
{
  part = ii[part.start]
  for(i in (n.days+1):1){
    out.trace[,part.start,i] <- out$state[,part,i]
    part <- out$parent[part,i]
  }
  print(part.start)
}
save(out.trace, file = paste0("./output/covid_pf_API_trace_", region, "_", n, "-", Sys.Date(), ".rdata"))
