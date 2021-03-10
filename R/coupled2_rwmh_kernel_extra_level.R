#' @rdname coupled2_rwmh_kernel_extra_level
#' @title RWM kernel for two coupled chains that target distribitions at different levels
#' @description generation of coupled proposals + accept/reject step
#' @param level an integer that determines probability distribution
#' @param state1 a list with state of the first chain
#' @param state2 a list with state of the second chain
#' @param identical a boolean variable that is True if chains are identical and False otherwise (it is always False for this function and it is left only for compatibility with other couplings)
#' @param tuning a list of parameters neeeded for RWM (standard deviation)
#' @param proposal_coupling function that determing proposal coupling
#'@return a list that contains states of two chains, updated value of the flag "identical" and cost of computations
#'@export

coupled2_rwmh_kernel_extra_level <- function( level,
                                              state1, state2,
                                              identical,
                                              tuning,
                                              proposal_coupling){
  cost = 0  # numerical cost of calculations
  # extract states
  chain_state1 <- state1$chain_state
  chain_state2 <- state2$chain_state

  #print("chain_state1")
  #print(chain_state1)
  #print("chain_state2")
  #print(chain_state2)
  # extract pdfs
  current_pdf1 <- state1$current_pdf
  current_pdf2 <- state2$current_pdf

  # extract tuning parameters
  proposal_sd <- tuning$proposal_sd

  # propose from 2-way coupling (output list with state1, state2, identical)
  proposal_value <- proposal_coupling(chain_state1, chain_state2,
                                      identical, tuning)
  cost = cost + proposal_value$cost
  # evaluate target density at proposals on same level
  proposal_state1 <- proposal_value$state1
  proposal_pdf1 <- logtarget(level - 1, proposal_state1)
  cost = cost + 2 ^ (level - 1)
  proposal_state2 <- proposal_value$state2
  proposal_pdf2 <- logtarget(level, proposal_state2)
  cost = cost + 2 ^ level
  

  # accept or reject proposals
  logu <- log(runif(1))
  accept1 <- FALSE
  accept2 <- FALSE

  if (is.finite(proposal_pdf1)) accept1 <- (logu < (proposal_pdf1 - current_pdf1))
  if (is.finite(proposal_pdf2)) accept2 <- (logu < (proposal_pdf2 - current_pdf2))
  #print("accept1 value")
  #print(accept1)
  #print("accept2 value")
  #print(accept2)

  if(accept1){
    chain_state1 <- proposal_state1
    current_pdf1 <- proposal_pdf1
  }
  if(accept2){
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

