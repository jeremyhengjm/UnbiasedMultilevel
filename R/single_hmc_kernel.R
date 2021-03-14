#' @rdname single_kernel_hmc
#' @title HMC kernel for a single MCMC chain that target probability distribution at the given level
#' @description Generation of HMC proposal vil leapfrog integration and accept/reject step
#' @param level a integer that determines density of probability distribution inthe multi-level approach
#' @param state a list with current position of the particle (element of the chain) and log of the density of probability distribution
#' @param tuning a list of parameters for HMC iteration: stepsize and number of steps
#'@return a updated state of the chain and cost of computations
#'@export

single_hmc_kernel <- function(level, state, tuning){

  # extract state and pdf
  chain_state <- state$chain_state
  current_pdf <- state$current_pdf

  # sample velocity or momentum
  current_v <- rnorm(dimension)

  # run leapfrog integrator
  leapfrog_result <- leapfrog(level, chain_state, current_v, tuning)
  proposed_v <- - leapfrog_result$v
  proposed_x <- leapfrog_result$x
  cost <- leapfrog_result$cost

  # compute pdf of proposal
  proposed_pdf <- logtarget(level, proposed_x)
  cost <- cost + 2^level

  # compute acceptance probability
  logacceptprob <- proposed_pdf - current_pdf
  logacceptprob <- logacceptprob + sum(current_v^2) / 2 - sum(proposed_v^2) / 2

  # accept or reject proposal
  accept <- FALSE
  if (is.finite(logacceptprob)){
    accept <- (log(runif(1)) < logacceptprob)
  }

  if (accept){
    chain_state <- proposed_x
    current_pdf <- proposed_pdf
  }

  chain_data = list(chain_state = chain_state, current_pdf = current_pdf, accept = accept)
  return(list(chain_data = chain_data, cost = cost))

}

