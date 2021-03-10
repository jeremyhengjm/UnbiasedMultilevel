#' @rdname coupled2_hmc_kernel_extra_level
#' @title HMC kernel for two coupled chains that target distribitions at different levels
#' @description generation of coupled proposals + accept/reject step
#' @param level an integer that determines probability distribution
#' @param state1 a list with state of the first chain
#' @param state2 a list with state of the second chain
#' @param identical a boolean variable that is True if chains are identical and False otherwise (it is always False for this function and it is left only for compatibility with other couplings)
#' @param tuning a list of parameters neeeded for HMC (aka time step for leapfrog integration)
#' @param proposal_coupling function that determing proposal coupling
#'@return a list that contains states of two chains, updated value of the flag "identical" and cost of computations
#'@export

coupled2_hmc_kernel_extra_level <- function(level, state1, state2, identical, tuning, proposal_coupling){

  cost = 0  # n umerical cost of calculations
  # extract state and pdf
  chain_state1 <- state1$chain_state
  chain_state2 <- state2$chain_state

  current_pdf1 <- state1$current_pdf
  current_pdf2 <- state2$current_pdf

  # propose from 2-way coupling
  # (output list with state1, velocity1, initial_velocity1,
  # state2, velocity2, initial_velocity2, identical)
  proposal_value <- proposal_coupling(level, chain_state1, chain_state2, identical, tuning)

  # evaluate target density at proposals on same level
  proposal_state1 <- proposal_value$state1
  cost = cost + proposal_value$cost
  proposal_velocity1 <- proposal_value$velocity1
  initial_velocity1 <- proposal_value$initial_velocity1
  proposal_pdf1 <- logtarget(level - 1, proposal_state1)
  cost = cost + 2 ^ (level - 1)

  proposal_state2 <- proposal_value$state2
  proposal_velocity2 <- proposal_value$velocity2
  initial_velocity2 <- proposal_value$initial_velocity2
  proposal_pdf2 <- logtarget(level, proposal_state2)
  cost = cost + 2 ^ level


  # compute acceptance probability
  logacceptprob1 <- proposal_pdf1 - current_pdf1
  logacceptprob1 <- logacceptprob1 + sum(initial_velocity1^2) / 2 - sum(proposal_velocity1^2) / 2

  logacceptprob2 <- proposal_pdf2 - current_pdf2
  logacceptprob2 <- logacceptprob2 + sum(initial_velocity2^2) / 2 - sum(proposal_velocity2^2) / 2

  # accept or reject proposals
  logu <- log(runif(1))
  accept1 <- FALSE
  accept2 <- FALSE
  if (is.finite(logacceptprob1)){
    accept1 <- (logu < logacceptprob1)
  }
  if (is.finite(logacceptprob2)){
    accept2 <- (logu < logacceptprob2)
  }

  if (accept1){
    chain_state1 <- proposal_state1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal_state2
    current_pdf2 <- proposal_pdf2
  }

  # check if chains are identical
  if (!identical){
    identical <- (proposal_value$identical && accept1 && accept2)
  }

  return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
              state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
              identical = identical, cost = cost))
}

