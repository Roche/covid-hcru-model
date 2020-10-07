#' A function to perform the kernel density particle filter
#'
#' @param y a matrix with no, dimension of each observation, rows and  nt, number of observation time points, columns
#' @param dllik a function to evaluate the logarithm of the likelihood with arguments y, x, and theta for data, state, and fixed parameters respectively
#' @param pstate a function to produce a point estimate of the state at time t+1 given the state at time t and fixed parameters theta, often the expectation
#' @param revo a function to propagate the state with arguments x, the current state, and theta, the fixed parameters
#' @param rprior a function to sample from the prior for the state and fixed parameters returns a list with elements x and theta
#' @param delta the kernel density smoothing parameter in (0,1) usually, around 0.99
#' @param n the number of particles
#' @param progress a boolean to display progress bar if TRUE
#' @param ... arguments passed on to resample
#' @return a list containing an n x (nt+1) matrix of normalized particle weights, a np x n x (nt+1) array of theta draws, a ns x n x (nt+1) array of state draws, and an n x nt parent matrix
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @export kd_pf
#' @references Liu, J. and West, M. Combined parameter and state estimation in 
#'   simulation-based filtering. in Sequential Monte Carlo Methods in Practice. 
#'   Doucet, A.; De Freitas, J. F. G. & Gordon, N. J. (Eds.), Springer-Verlag, 
#'   2001, 197-217
#' @seealso \code{\link{resample}}
#'
kd_pf = function(y, dllik, pstate, revo, rprior, n, delta=0.99, progress = TRUE, print = FALSE, ...)
{
  # Require R package to do resampling
  require(smcUtils)
  
  # Get number of time points
  if (!is.matrix(y)) y = matrix(y, 1)
  no = nrow(y) # not currently used
  nt = ncol(y)
  
  # Find dimension of state
  current.seed = .Random.seed
  tmp = rprior()
  ns = length(tmp$x)
  np = length(tmp$theta)
  .Random.seed = current.seed
  
  # Set up initial state
  state = array(NA, dim=c(ns,n,nt+1))
  theta = array(NA, dim=c(np,n,nt+1))
  for (j in 1:n) 
  {
    tmp = rprior()
    state[,j,1] = tmp$x
    theta[,j,1] = tmp$theta
  }
  
  # Initialize weights
  weight = matrix(NA, n, nt+1)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)
  
  # Kernel density smoothing parameter
  a = (3*delta-1)/(2*delta)
  h = sqrt(1-a^2)
  
  # Run particle filter
  if(progress & !print) pb = txtProgressBar(0,nt,style=3)
  p.state = matrix(NA, ns, n)
  p.theta = matrix(NA, np, n)
  p.weights = numeric(n)
  for (i in 1:nt) 
  {
    if(progress & !print) setTxtProgressBar(pb,i)

    # Kernel density estimate
    if(dim(theta)[1] == 1) m.theta = as.matrix(theta[,,i]) else m.theta = t(theta[,,i])
    theta.est = cov.wt(m.theta, weight[,i])
    mn = theta.est$center
    vr = theta.est$cov
    cl = try(t(chol(vr)),silent=TRUE) # If particles have degenerated on fixed parameters, 
    if(class(cl) == "try-error") cl <- suppressWarnings(t(chol(vr, pivot=TRUE))) # Cholesky decomp. may not be positive definite
    for (j in 1:n)
    {
      p.theta[,j] = a*theta[,j,i] +(1-a)*mn
    }
    
    # Look ahead
    for (j in 1:n)
    {
      p.state[,j]  = pstate(state[,j,i], theta[,j,i])
      p.weights[j] = log(weight[j,i]) + dllik(y[,i], p.state[,j], p.theta[,j])
    }
    
    # Normalize and resample auxiliary weights
    p.weights = renormalize(p.weights, log=T)
    tmp = resample(p.weights, log=FALSE, ...)
    kk = tmp$indices
    did.resample = !(all(kk == 1:n))

    # Regenerate fixed parameters, propagate state, and calculate weights
    for (j in 1:n) 
    {
      if(did.resample) theta[,j,i+1] = p.theta[,kk[j]] + h*cl%*%rnorm(np) else theta[,j,i+1] = theta[,j,i]
      state[,j,i+1] = revo(state[,kk[j],i], theta[,j,i+1])
      weight[j,i+1] = log(tmp$weights[kk[j]]) +
        dllik(y[,i],   state[,   j,i+1], theta[,   j,i+1]) -
        dllik(y[,i], p.state[,kk[j]],  p.theta[,kk[j]])
    }
    
    # Normalize weights
    weight[,i+1] = renormalize(weight[,i+1], log=T)
    if(!(all(weight[,i+1] < 1))) stop("Degeneration to a single particle")
    
    # Save parent for particle trace
    parent[,i+1] = kk
    
    if(print) print(i)
  }
  
  return(list(state=state, theta=theta, weight=weight, parent=parent))
}
