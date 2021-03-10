#' @rdname coupled2_rwmh_kernel_ext
#' @title RWM kernel for two chains that target distributions at different levels with synchronous couling of pCN proposals
#' @description generation of coupled proposals + accept/reject step
#' @param level an integer that determines probability distribution
#' @param state1 a list with state of the first chain
#' @param state2 a list with state of the second chain
#' @param tuning a list of parameters neeeded for RWM (standatd deviation)
#'@return a list that contains states of two chains and cost of computations
#'@export

coupled2_rwmh_kernel_ext <- function(level, state1, state2, tuning1, tuning2){
  cost = 0  # numerical cost of calculations
  # extract state and pdf
  chain_state1 <- state1$chain_state
  chain_state2 <- state2$chain_state

  current_pdf1 <- state1$current_pdf
  current_pdf2 <- state2$current_pdf

  # propose from 2-way coupling
  # (output list with state1, velocity1, initial_velocity1,
  # state2, velocity2, initial_velocity2, identical)
  proposal_value <- synchronous2_rwmh_coupling_ext(level, chain_state1, chain_state2, tuning1, tuning2)
  cost = cost + proposal_value$cost

  # evaluate target density at proposals on same level
  proposal_state1 <- proposal_value$state1
  proposal_pdf1 <- logtarget(level - 1, proposal_state1)
  cost = cost + 2 ^ (level - 1)
  
  proposal_state2 <- proposal_value$state2
  proposal_pdf2 <- logtarget(level, proposal_state2)
  cost = cost + 2 ^ level
  
  # compute acceptance probability
  logacceptprob1 <- proposal_pdf1 - current_pdf1
  logacceptprob2 <- proposal_pdf2 - current_pdf2

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
  
  return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
              state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
              cost = cost))
}


