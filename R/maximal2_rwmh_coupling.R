#' @rdname maximal2_rwmh_coupling
#' @title maximal coupling of 2 RWM proposals
#' @description Generate pCn proposals for two chains via maximal coupling
#' @param chain_state1 a vector with coordinates of the first particle
#' @param chain_state2 a vector with coordinates of the second particle
#' @param identical a flag that is True if chains are identical and False otherwise
#' @param tuning a list that contains parameters needed for RWM: standard devation
#'@return a list that contains state of the first chain, state of the second chain, updated value of the flag "identical", cost of proposal generation
#'@export

maximal2_rwmh_coupling <- function(chain_state1, chain_state2, identical, tuning){
  cost <- 0  # cost of proposal generation
  # extract tuning parameters
  proposal_sd <- tuning$proposal_sd

  # sample first proposal
  state1 <- chain_state1 + proposal_sd * rnorm(dimension)

  # evaluate proposal transition densities at first proposal
  pdf1 <- sum(dnorm(state1, chain_state1, proposal_sd, log = TRUE))
  pdf2 <- sum(dnorm(state1, chain_state2, proposal_sd, log = TRUE))
  logacceptprob <- min(pdf1, pdf2) - pdf1

  if (log(runif(1)) < logacceptprob){
    # return common proposal for both chains
    return(list(state1 = state1, state2 = state1, identical = TRUE, cost = cost))

  } else {
    # rejection sampler for second proposal
    reject <- TRUE
    state2 <- NA
    while (reject){
      # sample proposed value
      state2 <- chain_state2 + proposal_sd * rnorm(dimension)

      # evaluate proposal transition densities at proposed value
      pdf1 <- sum(dnorm(state2, chain_state1, proposal_sd, log = TRUE))
      pdf2 <- sum(dnorm(state2, chain_state2, proposal_sd, log = TRUE))
      logrejectprob <- min(pdf1, pdf2) - pdf2

      # accept or reject
      reject <- (log(runif(1)) < logrejectprob)

    }

    return(list(state1 = state1, state2 = state2, identical = FALSE, cost = cost))

  }

}
