# covid_functions.r - functions needed for running particle filters for the COVID-19 pandemic model based off Moghadas et al. 2020

# covid_prop - function that propagates the state of the epidemic from one time point to the next
# disease states: S, E, Aasym, Aasymiso, Amild, Amildiso, Isev, Iseviso, Isevhosp, Isevhospiso, Tasym, Tmild, Tsev, Tsevhosp, Hfloor, Hicu, Hvent, D, R, Rdiff, Inc, and inc are lambda-length vectors
# social-distancing / change in testing over time: xi, phi_x, sigma_w, w_A, phi_A, w_I, phi_I are scalars
# population size: P is a lambda-length vector
# contact matrices: M, Mtilde are lambda by lambda matrices
# age-dependent probabilities: zeta, theta, rho, c are lambda-length vectors
# transmission parameters: beta, kappa_M, kappa_A, sigma, gamma, f_A, f_I, tau_A, tau_I, q are scalars
# testing parameters: w_s, r_A, r_I, f_T, phi, omega are scalars
# hospital parameters: delta, v, m_h, m_c, m_v, mu_h, mu_c, mu_v, psi_h, psi_c, psi_v are scalars
covid_prop = function(S, E, Aasym, Aasymiso, Amild, Amildiso, Isev, Iseviso, Isevhosp, Isevhospiso, # states of the model - transmission
                      Tasym, Tmild, Tsev, Tsevhosp, Hfloor, Hicu, Hvent, D, R, # states of the model - resource use
                      xi, # state of the model - local adjustment to transmission / social distancing
                      phi_x, sigma_x, # AR(1) parameters for social distancing
                      w_A, w_I, # states of the model - proportion of symptomatic individuals that get tested
                      phi_A, phi_I, sigma_w, # AR(1) parameters for proportion tested
                      M, Mtilde, P, # fixed constants
                      zeta, theta, rho, c, # age-dependent parameter vectors
                      beta, kappa_M, kappa_A, sigma, gamma, # transmission parameters
                      f_A, f_I, tau_A, tau_I, q, # behavioral parameters
                      w_s, r_A, r_I, f_T, phi, omega, alpha, # testing parameters
                      delta, v, m_h, m_c, m_v, mu_h, mu_c, mu_v, psi_h, psi_c, psi_v,  # resource-use parameters
                      random = TRUE, xi.fixed = NULL, w_A.fixed = NULL, w_I.fixed = NULL)
{
  # Check argument dimensions
  lambda <- length(S)
  stopifnot(length(Isev) == lambda & length(Isevhosp) == lambda & length(Amild) == lambda & length(Aasym) == lambda & length(P) == lambda)
  stopifnot(length(Tasym) == lambda & length(Tmild) == lambda & length(Tsev) == lambda & length(Tsevhosp) == lambda)
  stopifnot(length(Hfloor) == lambda & length(Hicu) == lambda & length(Hvent) == lambda)
  stopifnot(length(D) == lambda & length(R) == lambda)
  stopifnot(length(Iseviso) == lambda & length(Isevhospiso) == lambda & length(Aasymiso) == lambda & length(Amildiso) == lambda)
  stopifnot(all(dim(M) == lambda) & all(dim(Mtilde) == lambda) & length(dim(M)) == 2 & length(dim(Mtilde)) == 2)
  stopifnot(length(zeta) == lambda & length(theta) == lambda & length(rho) == lambda & length(c) == lambda)
  stopifnot(length(beta) == 1 & length(alpha) == 1 & length(kappa_M) == 1 & length(kappa_A) == 1 & length(sigma) == 1 & length(gamma) == 1)
  stopifnot(length(f_A) == 1 & length(f_I) == 1 & length(tau_A) == 1 & length(tau_I) == 1 & length(q) == 1)
  stopifnot(length(w_s) == 1 & length(r_A) == 1 & length(r_I) == 1 & length(f_T) == 1 & length(phi) == 1 & length(omega) == 1)
  stopifnot(length(v) == 1 & length(m_h) == 1 & length(m_c) == 1 & length(m_v) == 1 & length(mu_h) == 1 & length(mu_c) == 1 & length(mu_v) == 1)
  stopifnot(length(psi_h) == 1 & length(psi_c) == 1 & length(psi_v) == 1 & length(delta) == 1)
  stopifnot(length(xi) == 1 & length(phi_x) == 1 & length(sigma_x) == 1)
  stopifnot(length(w_I) == 1 & length(phi_A) == 1 & length(sigma_w) == 1)
  stopifnot(length(w_A) == 1 & length(phi_I) == 1)
  stopifnot(is.null(xi.fixed) | xi %in% 1:length(xi.fixed))
  stopifnot(is.null(w_A.fixed) | w_A %in% 1:length(w_A.fixed))
  stopifnot(is.null(w_I.fixed) | w_I %in% 1:length(w_I.fixed))

  # Incidence vector operations
  if(is.null(xi.fixed)){
    t_rate <- beta*xi*S 
  } else {
    t_rate <- beta*xi.fixed[xi]*S
  }
  
  # Incidence matrix operations
  p_comm <- (Isev + Isevhosp + kappa_M*Amild + kappa_A*Aasym + (1 - f_T)*(Tsev + Tsevhosp + kappa_M*Tmild + kappa_A*Tasym)) / P
  comm_contacts <- M %*% p_comm
  p_iso <- (Iseviso + Isevhospiso + kappa_M*Amildiso + kappa_A*Aasymiso + f_T*(Tsev + Tsevhosp + kappa_M*Tmild + kappa_A*Tasym)) / P
  iso_contacts <- Mtilde %*% p_iso
  inc <- t_rate*(comm_contacts[,1] + iso_contacts[,1])
  
  # Propagate susceptibles and exposed
  S1 <- S - inc
  E1 <- E + inc - sigma*E
  
  # Set proportion of testing throughput limits for each testing / age category
  if(sum(Tasym, Tmild, Tsev, Tsevhosp) <= 0){
    asym_fac <- mild_fac <- sev_fac <- hosp_fac <- 0
  } else {
    asym_fac <- Tasym / sum(Tasym, Tmild, Tsev, Tsevhosp)
    mild_fac <- Tmild / sum(Tasym, Tmild, Tsev, Tsevhosp)
    sev_fac <- Tsev / sum(Tasym, Tmild, Tsev, Tsevhosp)
    hosp_fac <- Tsevhosp / sum(Tasym, Tmild, Tsev, Tsevhosp)
  }
  
  # Set proportion tested and testing limit parameters to fixed trace if it is provided
  if(!is.null(w_A.fixed)){
    w_At <- w_A
    w_A <- w_A.fixed[w_At]
  }
  if(!is.null(w_I.fixed)){
    w_It <- w_I
    w_I <- w_I.fixed[w_It]
  }

  # Propagate asymptomatic individuals
  Aasym1 <- Aasym + zeta*sigma*E - (1 - w_s)*gamma*Aasym - w_s*r_A*Aasym + 
    phi*(1-omega)*apply(rbind(Tasym, alpha * asym_fac), 2, min)
  Aasymiso1 <- Aasymiso - gamma*Aasymiso + phi*omega*apply(rbind(Tasym, alpha * asym_fac), 2, min)
  
  # Propagate mildly symptomatic individuals
  Amild1 <- Amild + (1 - zeta)*theta*sigma*E - (1 - w_A)*(1 - f_A)*gamma*Amild - (1 - w_A)*f_A*tau_A*Amild -
    w_A*r_A*Amild + phi*(1-omega)*apply(rbind(Tmild, alpha * mild_fac), 2, min)
  Amildiso1 <- Amildiso + (1 - w_A)*f_A*tau_A*Amild - gamma*Amildiso + 
    phi*omega*apply(rbind(Tmild, alpha * mild_fac), 2, min)
  
  # Propagate severe symptomatic individuals
  Isev1 <- Isev + (1 - zeta)*(1 - theta)*(1 - q)*(1 - rho)*sigma*E - (1 - w_I)*(1 - f_I)*gamma*Isev - 
    (1 - w_I)*f_I*tau_I*Isev - w_I*r_I*Isev + 
    phi*(1-omega)*apply(rbind(Tsev, alpha * sev_fac), 2, min)
  Iseviso1 <- Iseviso + (1 - zeta)*(1 - theta)*q*(1 - rho)*sigma*E - gamma*Iseviso + (1 - w_I)*f_I*tau_I*Isev +
    phi*omega*apply(rbind(Tsev, alpha * sev_fac), 2, min)
  Isevhosp1 <- Isevhosp + (1 - zeta)*(1 - theta)*(1 - q)*rho*sigma*E - (1 - w_I)*(1 - f_I)*delta*Isevhosp - 
    (1 - w_I)*f_I*tau_I*Isevhosp - w_I*r_I*Isevhosp + 
    phi*(1-omega)*apply(rbind(Tsevhosp, alpha * hosp_fac), 2, min)
  Isevhospiso1 <- Isevhospiso + (1 - zeta)*(1 - theta)*q*rho*sigma*E - delta*Isevhospiso + (1 - w_I)*f_I*tau_I*Isevhosp + 
    phi*omega*apply(rbind(Tsevhosp, alpha * hosp_fac), 2, min)
  
  # Propagate individuals under testing
  Tasym1 <- Tasym + w_s*r_A*Aasym - phi*apply(rbind(Tasym, alpha * asym_fac), 2, min) - gamma*Tasym
  Tasympos <- phi*omega*apply(rbind(Tasym, alpha * asym_fac), 2, min)
  Tmild1 <- Tmild + w_A*r_A*Amild - phi*apply(rbind(Tmild, alpha * mild_fac), 2, min) - gamma*Tmild
  Tmildpos <- phi*omega*apply(rbind(Tmild, alpha * mild_fac), 2, min)
  Tsev1 <- Tsev + w_I*r_I*Isev - phi*apply(rbind(Tsev, alpha * sev_fac), 2, min) - gamma*Tsev
  Tsevpos <- phi*omega*apply(rbind(Tsev, alpha * sev_fac), 2, min)
  Tsevhosp1 <- Tsevhosp + w_I*r_I*Isevhosp - phi*apply(rbind(Tsevhosp, alpha * hosp_fac), 2, min) - delta*Tsevhosp
  Thosppos <- phi*omega*apply(rbind(Tsevhosp, alpha * hosp_fac), 2, min)
  
  # Calculate new hospitalizations
  newfloor <- (1-c)*(1-w_I)*(1-f_I)*delta*Isevhosp + (1-c)*delta*Isevhospiso + (1-c)*delta*Tsevhosp
  newicu <- c*(1-v)*(1-w_I)*(1-f_I)*delta*Isevhosp + c*(1-v)*delta*Isevhospiso + c*(1-v)*delta*Tsevhosp
  newvent <- c*v*(1-w_I)*(1-f_I)*delta*Isevhosp + c*v*delta*Isevhospiso + c*v*delta*Tsevhosp
  
  # Calculate new deaths
  deathfloor <- m_h*mu_h*Hfloor
  deathicu <- m_c*mu_c*Hicu
  deathvent <- m_v*mu_v*Hvent
  
  # Calculate new recoveries from hospital
  recfloor <- (1-m_h)*psi_h*Hfloor
  recicu <- (1-m_c)*psi_c*Hicu
  recvent <- (1-m_v)*psi_v*Hvent
  
  # Propagate individuals in the hospital
  Hfloor1 <- Hfloor + newfloor - deathfloor - recfloor
  Hicu1 <- Hicu + newicu - deathicu - recicu
  Hvent1 <- Hvent + newvent - deathvent - recvent
  
  # Track cumulative incidence, incidence, deaths, and recoveries
  D1 <- D + deathfloor + deathicu + deathvent
  R1 <- R + recfloor + recicu + recvent +
    (1 - w_s)*gamma*Aasym + gamma*Aasymiso + 
    (1 - w_A)*(1 - f_A)*gamma*Amild + gamma*Amildiso +
    (1 - w_I)*(1 - f_I)*gamma*Isev + gamma*Iseviso +
    gamma*Tasym + gamma*Tmild + gamma*Tsev
  inc1 <- inc
  Tpos1 <- Tasympos + Tmildpos + Tsevpos + Thosppos
  Hnew1 <- newfloor + newicu + newvent
  
  # Propagate local adjustment and testing parameters
  if(is.null(xi.fixed)) xi1 <- max(phi_x*xi + random * rnorm(1, 0, sigma_x), 0) else xi1 <- xi + 1
  if(is.null(w_A.fixed) | is.null(w_I.fixed)){
    e <- rnorm(1, 0, sigma_w)
    if(is.null(w_A.fixed)) w_A1 <- min(max(phi_A*(w_A - 1) + random * e + 1, 0), 1) else w_A1 <- w_At + 1
    if(is.null(w_I.fixed)) w_I1 <- min(max(phi_I*(w_I - 1) + random * e + 1, 0), 1) else w_I1 <- w_It + 1
  } else {
    w_A1 <- w_At + 1
    w_I1 <- w_It + 1
  }

  # Track rounding errors  
  Rcheck <- P - (S1 + E1 + Isev1 + Isevhosp1 + Iseviso1 + Isevhospiso1 + Amild1 + Amildiso1 + Aasym1 + Aasymiso1 + Tasym1 + Tmild1 + Tsev1 + Tsevhosp1 + Hfloor1 + Hicu1 + Hvent1 + D1)
  
  # Check that the population is closed
  tot_curr <- S + E + Isev + Isevhosp + Iseviso + Isevhospiso + Amild + Amildiso + Aasym + Aasymiso + Tasym + Tmild + Tsev + Tsevhosp + Hfloor + Hicu + Hvent + D + R
  tot_next <- S1 + E1 + Isev1 + Isevhosp1 + Iseviso1 + Isevhospiso1 + Amild1 + Amildiso1 + Aasym1 + Aasymiso1 + Tasym1 + Tmild1 + Tsev1 + Tsevhosp1 + Hfloor1 + Hicu1 + Hvent1 + D1 + R1
  stopifnot(all(abs(P - tot_curr) < 1e-6) & all(abs(P - tot_next) < 1e-6))
    
  return(list(S=S1, E=E1, Aasym=Aasym1, Aasymiso=Aasymiso1, Amild=Amild1, Amildiso=Amildiso1, Isev=Isev1, Iseviso=Iseviso1, Isevhosp=Isevhosp1, Isevhospiso=Isevhospiso1,
              Tasym=Tasym1, Tmild=Tmild1, Tsev=Tsev1, Tsevhosp=Tsevhosp1, Hfloor=Hfloor1, Hicu=Hicu1, Hvent=Hvent1, D=D1, R=R1,
              inc=inc1, Tpos=Tpos1, Hnew = Hnew1,
              xi = xi1, w_A = w_A1, w_I = w_I1,
              Rdiff=R1-Rcheck))
}

# dllik - function to return the log of the likelihood of observed pandemic data at a given time given the current state of pandemic and fixed parameters
# Tpos, Hfloor, Hicu, Hvent, D are lambda-length vectors
# sigma_pos, sigma_hosp, sigma_icu, sigma_vent, sigma_d are scalars
# loglik determines which parts of the likelihood function are evaluated
covid_dllik = function(y, Tpos, Hfloor, Hicu, Hvent, D,
                 sigma_pos, sigma_hosp, sigma_icu, sigma_vent, sigma_d,
                 loglik = c('t', 'h', 'i', 'v', 'd'))
{
  # Calculate log-likelihood of observations and set to NA if zero observed
  if(sum(Tpos) <= 0 | !('t' %in% loglik)) tloglik <- NA else tloglik <- dlnorm(y[1], log(sum(Tpos)), sigma_pos, log = TRUE)
  if(sum(Hfloor + Hicu + Hvent) <= 0 | !('h' %in% loglik)) hloglik <- NA else hloglik <- dlnorm(y[2], log(sum(Hfloor + Hicu +  Hvent)), sigma_hosp, log=TRUE)
  if(sum(Hicu) <= 0 | !('i' %in% loglik)) iloglik <- NA else iloglik <- dlnorm(y[3], log(sum(Hicu + Hvent)), sigma_icu, log=TRUE)
  if(sum(Hvent) <= 0 | !('v' %in% loglik)) vloglik <- NA else vloglik <- dlnorm(y[4], log(sum(Hvent)), sigma_vent, log = TRUE)
  if(sum(D) <= 0 | !('d' %in% loglik)) dloglik <- NA else dloglik <- dlnorm(y[5], log(sum(D)), sigma_d, log = TRUE)

  # Debug
  if(any(c(tloglik, hloglik, iloglik, vloglik, dloglik) %in% c(-Inf, Inf))) stop("Infinite log-likelihood")

  # Compute log-likelihood
  return(sum(tloglik, hloglik, iloglik, vloglik, dloglik, na.rm = TRUE))
}

# rprior - function that samples from the prior distribution of the fixed parameters and the initial state of the pandemic
rprior = function(loglik = c('t', 'h', 'i', 'v', 'd'), P = 330000000 * c(.25, .4, .19, .16),
                  xi.fixed = NULL, w_A.fixed = NULL, w_I.fixed = NULL, p_a.fixed = NULL,
                  theta.post = NULL, weight = NULL)
{
  # Sample from priors for fixed parameters
  if(is.null(theta.post)){ # If no prior samples input, sample from prespecified prior distributions
    theta <- c()
    
    if('t' %in% loglik){
      p_alpha.mean <- 0.1 # proportion of testing capacity applied to infected population
      sigma_pos.mean <- 0.5 # standard deviation of the log of new positive tests
      theta <- c(theta,
                 p_alpha = u2theta(rbeta(1, 50 * p_alpha.mean, 50*(1 - p_alpha.mean)), 0, 1), # p_alpha
                 sigma_pos = log(rgamma(1, 50 * sigma_pos.mean, 50)) # sigma_pos
      )
    }
    
    if('h' %in% loglik){
      rho.mean <- c(0.592, 0.350, 0.573, 0.698) # Proportion of severe cases requiring hospitalization
      sigma_hosp.mean <- .5 # standard deviation of log of currently hospitalized
      delta.mean <- 1/7 # time to hospitalization from symptom onset
      psi_h.mean <- 1/8.3 # rate of recovery from hospital floor
      psi_c.mean <- 1/22 # rate of recovery from ICU w/o vent
      psi_v.mean <- 1/28 # rate of recovery from ICU w/ vent
      theta <- c(theta,
                 rho = u2theta(rbeta(length(P), 50*rho.mean, 50*(1 - rho.mean)), 0, 1), # rho
                 sigma_hosp = log(rgamma(1, 50 * sigma_hosp.mean, 50)), # sigma_h
                 delta = log(rgamma(1, 50 * delta.mean, 50)), # delta
                 psi_h = log(rgamma(1, 50 * psi_h.mean, 50)), # psi_h
                 psi_c = log(rgamma(1, 50 * psi_c.mean, 50)), # psi_c
                 psi_v = log(rgamma(1, 50 * psi_v.mean, 50)) # psi_v
      )
    }
    
    if('i' %in% loglik){
      c.mean <- c(0.152, 0.215, 0.215, 0.237) # Proportion of hospitalized cases requiring admission to ICU
      sigma_icu.mean <- .5 # standard deviation of log of currently in ICU
      theta <- c(theta,
                 c = u2theta(rbeta(length(P), 10*c.mean, 10*(1 - c.mean)), 0, 1), # c
                 sigma_icu = log(rgamma(1, 50 * sigma_icu.mean, 50)) # sigma_icu
      )
    }
    
    if('v' %in% loglik){
      v.mean <- 0.6 # proportion of ICU cases requiring mechanical ventilation
      sigma_vent.mean <- .5 # standard deviation of log of currently on ventilation
      theta <- c(theta,
                 v = u2theta(rbeta(1, 10*v.mean, 10*(1 - v.mean)), 0, 1), # v
                 sigma_vent = log(rgamma(1, 50 * sigma_vent.mean, 50)) # sigma_vent
      )
    }
    
    if('d' %in% loglik){
      m_h.mean <- 0.067 # weight associated with death on the hospital floor
      m_c.mean <- 0.048 # weight associated with death in the ICU without ventilation
      m_v.mean <- 0.099 # weight associated with death in the ICU for ventilated patients
      mu_h.mean <- 1 / 9.7 # inverse of the time until death on hospital floor
      mu_c.mean <- 1 / 6 # inverse of the time until death in the ICU (non-ventilated)
      mu_v.mean <- 1 / 5.4 # inverse of the time until death in the ICU (ventilated)
      sigma_d.mean <- 0.5 # standard deviation of log of observed deaths
      theta <- c(theta,
                 m_h = u2theta(rbeta(1, 10*m_h.mean, 10*(1 - m_h.mean)), 0, 1), # m_h
                 m_c = u2theta(rbeta(1, 10*m_c.mean, 10*(1 - m_c.mean)), 0, 1), # m_c
                 m_v = u2theta(rbeta(1, 10*m_v.mean, 10*(1 - m_v.mean)), 0, 1), # m_v
                 mu_h = log(rgamma(1, 50 * mu_h.mean, 50)), # mu_h
                 mu_c = log(rgamma(1, 50 * mu_c.mean, 50)), # mu_c
                 mu_v = log(rgamma(1, 50 * mu_v.mean, 50)), # mu_v
                 sigma_d = log(rgamma(1, 50 * sigma_d.mean, 50)) # sigma_d
      )
    }
  } else { # If prior samples provided, sample a particle randomly
    stopifnot(!is.null(weight)) # Need posterior sample weights
    n.part <- length(weight) # How many particles in posterior sample?
    stopifnot(n.part == dim(theta.post)[2]) # Number of particles should match length of weight vector
    
    # Sample a particle
    index <- sample(1:n.part, 1, prob = weight)
    theta <- theta.post[,index]
  }

  # Ensure finite values sampled
  stopifnot(all(!(theta %in% c(Inf, -Inf, NaN))))

  # Set initial state of pandemic
  lambda <- length(P)
  E0 <- rep(1, lambda)
  if(is.null(xi.fixed)) xi0 <- rlnorm(1, 2, .3) else xi0 <- 1 # Initial baseline social distancing / transmission adjustment
  if(is.null(w_A.fixed)) w_A0 <- rbeta(1, 100*.1, 100*(1-.1)) else w_A0 <- 1 # Initial baseline proportion of mild cases tested
  if(is.null(w_I.fixed)) w_I0 <- rbeta(1, 100*.1, 100*(1-.1)) else w_I0 <- 1 # Initial baseline proportion of severe cases tested
  x <- c(S = P - E0,
         E = E0, Isev = rep(0, lambda), Isevhosp = rep(0, lambda), Iseviso = rep(0, lambda), Isevhospiso = rep(0, lambda),
         Aasym = rep(0, lambda), Aasymiso = rep(0, lambda), Amild= rep(0, lambda), Amildiso = rep(0, lambda),
         Tasym = rep(0, lambda), Tmild = rep(0, lambda), Tsev = rep(0, lambda), Tsevhosp = rep(0, lambda), Tpos = rep(0, lambda),
         Hfloor = rep(0, lambda), Hicu = rep(0, lambda), Hvent = rep(0, lambda), Hnew = rep(0, lambda), D = rep(0, lambda), R = rep(0, lambda),
         xi = xi0, w_A = w_A0, w_I = w_I0)

  if(length(theta) > 0) return(list(theta=theta, x=x)) else return(x)
}

## Utility functions ##

# Functions to reparameterize theta to [a,b] from the real line and vice versa
theta2u = function(theta,a,b)
{
  etheta = exp(theta)
  return((b*etheta + a) / (1 + etheta))
}
u2theta = function(u,a,b)
{
  U = (u - a) / (b - a)
  return(log(U / (1 - U)))
}