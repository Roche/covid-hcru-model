#' A function to perform the auxiliary particle filter
#'
#' @param y a matrix with no, dimension of each observation, rows and  nt, number of observation time points, columns
#' @param dllik a function to evaluate the logarithm of the likelihood with arguments y and x for data and state respectively
#' @param pstate a function to produce a point estimate of the state at time t+1 given the state at time t, often the expectation
#' @param revo a function to propagate the state with argument the state x
#' @param rprior a function to sample for the prior for the state
#' @param n the number of particles
#' @param progress a boolean to display progress bar if TRUE
#' @param ... arguments passed on to resample
#' @return a list containing an n x (nt+1) matrix of normalized particles weights, a ns x n x (nt+1) array of state draws, and an n x nt parent matrix
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @export apf
#' @references Pitt, M.K.; Shephard, N. (1999). 
#'   "Filtering Via Simulation: Auxiliary Particle Filters". 
#'   Journal of the American Statistical Association. 94 (446): 590?591.
#' @seealso \code{\link{resample}}
#'
apf = function(y, dllik, pstate, revo, rprior, n, progress = TRUE, print = FALSE, ...)
{
  # Load R package to do resampling
  require(smcUtils)

  # Get number of time points
  if (!is.matrix(y)) y = matrix(y, 1)
  no = nrow(y) # not currently used
  nt = ncol(y)

  # Find dimension of state
  current.seed = .Random.seed
  ns = length(rprior())
  .Random.seed = current.seed

  # Set up initial state
  state = array(NA, dim=c(ns,n,nt+1))
  for (j in 1:n)
  {
    state[,j,1] = rprior()
  }

  # Initialize weights
  weight = matrix(NA, n, nt+1)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)

  # Run auxiliary particle filter
  if(progress & !print) pb = txtProgressBar(0,nt,style=3)
  p.state = matrix(NA, ns, n)
  p.weights = numeric(n)
  for (i in 1:nt) 
  {
    if(progress & !print) setTxtProgressBar(pb,i)

    # Look ahead
    for (j in 1:n)
    {
      p.state[,j] = pstate(state[,j,i])
      p.weights[j] = log(weight[j,i]) + dllik(y[,i], p.state[,j])
    }
    
    # Normalize and resample auxiliary weights
    p.weights = renormalize(p.weights, log=T)
    tmp = resample(p.weights, log=FALSE, ...)    
    kk = tmp$indices

    # Propagate state and calculate weights
    for (j in 1:n) 
    {
      state[,j,i+1] = revo(state[,kk[j],i])
      weight[j,i+1] = log(tmp$weights[kk[j]]) + 
                      dllik(y[,i],   state[,   j,i+1]) -
                      dllik(y[,i], p.state[,kk[j]])
    }

    # Normalize weights
    weight[,i+1] = renormalize(weight[,i+1], log=T)

    # Save parent for particle trace
    parent[,i+1] = kk
    
    if(print) print(i)
  }

  return(list(state=state, weight=weight, parent=parent))
}