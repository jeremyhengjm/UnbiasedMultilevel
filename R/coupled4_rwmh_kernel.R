#' @rdname coupled4_rwmh_kernel
#' @title RWM kernel for four coupled chains
#' @description generation of coupled proposals + accept/reject step
#' @param level an integer that determines probability distribution: level == l
#' @param state_coarse1 a list with state of the first chain at the level "l - 1"
#' @param state_coarse2 a list with state of the second chain at the level "l - 1"
#' @param state_fine1 a list with state of the first chain at the level "l"
#' @param state_fine2 a list with state of the second chain at the level "l"
#' @param identical_coarse a boolean variable that is True if chains at the level "l - 1" are identical and False otherwise
#' @param identical_fine a boolean variable that is True if chains at the level "l" are identical and False otherwise
#' @param tuning_coarse a list of parameters neeeded for proposal generation with RWM (aka standard deviation) for chains at the level "l - 1"
#' @param tuning_fine a list of parameters neeeded for proposal generation with RWM (aka standard deviation) for chains at the level "l"
#' @param proposal_coupling function that determing proposal coupling
#'@return a list that contains states all four chains, updated value of the flag "identical_coarse" and "identical_coarse", cost of computations
#'@export

coupled4_rwmh_kernel <- function(level,
                                 state_coarse1, state_coarse2,
                                 state_fine1, state_fine2,
                                 identical_coarse, identical_fine,
                                 tuning_coarse, tuning_fine,
                                 proposal_coupling){
  cost = 0  # numerical cost of calculations
  # extract states
  chain_state_coarse1 <- state_coarse1$chain_state
  chain_state_coarse2 <- state_coarse2$chain_state
  
  chain_state_fine1 <- state_fine1$chain_state
  chain_state_fine2 <- state_fine2$chain_state
  # print("the values has been copied")

  # extract pdfs
  current_pdf_coarse1 <- state_coarse1$current_pdf
  current_pdf_coarse2 <- state_coarse2$current_pdf

  current_pdf_fine1 <- state_fine1$current_pdf
  current_pdf_fine2 <- state_fine2$current_pdf

  # extract tuning parameters
  proposal_sd_coarse <- tuning_coarse$proposal_sd
  proposal_sd_fine <- tuning_fine$proposal_sd

  # propose from 4-way coupling (output list with coarse1, coarse2, fine1, fine2,
  # identical_coarse, identical_fine)
  proposal_value <- proposal_coupling(chain_state_coarse1, chain_state_coarse2,
                                      chain_state_fine1, chain_state_fine2,
                                      identical_coarse, identical_fine,
                                      tuning_coarse, tuning_fine)
  cost = cost + proposal_value$cost
  proposal_state_coarse1 <- proposal_value$coarse1
  proposal_pdf_coarse1 <- logtarget(level-1, proposal_state_coarse1)
  cost = cost + 2 ^ (level - 1)
  if (proposal_value$identical_coarse){
    proposal_state_coarse2 <- proposal_state_coarse1
    proposal_pdf_coarse2 <- proposal_pdf_coarse1
  } else {
    proposal_state_coarse2 <- proposal_value$coarse2
    proposal_pdf_coarse2 <- logtarget(level-1, proposal_state_coarse2)
    cost = cost + 2 ^ (level - 1)
  }

  # evaluate target density at proposals on fine level
  proposal_state_fine1 <- proposal_value$fine1
  proposal_pdf_fine1 <- logtarget(level, proposal_state_fine1)
  cost = cost + 2 ^ level
  if (proposal_value$identical_fine){
    proposal_state_fine2 <- proposal_state_fine1
    proposal_pdf_fine2 <- proposal_pdf_fine1
  } else {
    proposal_state_fine2 <- proposal_value$fine2
    proposal_pdf_fine2 <- logtarget(level, proposal_state_fine2)
    cost = cost + 2 ^ level
  }

  # accept or reject proposals
  logu <- log(runif(1))
  accept_coarse1 <- FALSE
  accept_coarse2 <- FALSE
  accept_fine1 <- FALSE
  accept_fine2 <- FALSE

  if (is.finite(proposal_pdf_coarse1)) accept_coarse1 <- (logu < (proposal_pdf_coarse1 - current_pdf_coarse1))
  if (is.finite(proposal_pdf_coarse2)) accept_coarse2 <- (logu < (proposal_pdf_coarse2 - current_pdf_coarse2))
  if (is.finite(proposal_pdf_fine1)) accept_fine1 <- (logu < (proposal_pdf_fine1 - current_pdf_fine1))
  if (is.finite(proposal_pdf_fine2)) accept_fine2 <- (logu < (proposal_pdf_fine2 - current_pdf_fine2))

  if (accept_coarse1){
    chain_state_coarse1 <- proposal_state_coarse1
    current_pdf_coarse1 <- proposal_pdf_coarse1
  }
  if (accept_coarse2){
    chain_state_coarse2 <- proposal_state_coarse2
    current_pdf_coarse2 <- proposal_pdf_coarse2
  }
  if (accept_fine1){
    chain_state_fine1 <- proposal_state_fine1
    current_pdf_fine1 <- proposal_pdf_fine1
  }
  if (accept_fine2){
    chain_state_fine2 <- proposal_state_fine2
    current_pdf_fine2 <- proposal_pdf_fine2
  }

  # check if chains are identical
  if (!identical_coarse){
    identical_coarse <- (proposal_value$identical_coarse && accept_coarse1 && accept_coarse2)
  }

  if (!identical_fine){
    identical_fine <- (proposal_value$identical_fine && accept_fine1 && accept_fine2)
  }

  return(list(state_coarse1 = list(chain_state = chain_state_coarse1, current_pdf = current_pdf_coarse1),
              state_coarse2 = list(chain_state = chain_state_coarse2, current_pdf = current_pdf_coarse2),
              state_fine1 = list(chain_state = chain_state_fine1, current_pdf = current_pdf_fine1),
              state_fine2 = list(chain_state = chain_state_fine2, current_pdf = current_pdf_fine2),
              identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))

}
