#' @rdname single_rwmh_kernel
#' @title RWM kernel for a single MCMC chain that target probability distribution at the given level
#' @description Generation of RWM proposal and accept/reject step
#' @param level a integer that determines density of probability distribution inthe multi-level approach
#' @param state a list with current position of the particle (element of the chain) and log of the density of probability distribution
#' @param tuning a list of parameters for RWM iteration: snadard devation
#'@return a updated state of the chain and cost of computations
#'@export

single_rwmh_kernel <- function(level, state, tuning){
  chain_state <- state$chain_state
  current_pdf <- state$current_pdf
  proposal_sd <- tuning$proposal_sd
  proposal_value <- proposal_sd * rnorm(dimension)
  #proposal_value <- 0.1 * rnorm(2)
  #print("dimension")
  #print(dimension)
  #print("proposal_value")
  #print(proposal_value[1])
  #print("proposal_sd")
  #print(proposal_sd)
  #print('tuning')
  #print(tuning)
  proposal_pdf <- logtarget(level, proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    chain_data = list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept)    
  } else {
    chain_data = list(chain_state = chain_state, current_pdf = current_pdf, accept = accept)
  }
  return(list(chain_data = chain_data, cost = 2 ^ level))
}

