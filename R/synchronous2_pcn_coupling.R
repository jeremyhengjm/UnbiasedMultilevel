#' @rdname synchronous2_pcn_coupling
#' @title Synchronous coupling of pCN proposals for two chains
#' @description Generation of two pCN proposals based on the common Gaussian increment
#' @param chain_state1 a vector with position of the first particle
#' @param chain_state2 a vector with position of the second particle
#' @param identical a boolean variable, that takes value "True" if particles are at the same position and 
#' @param tuning a list that contains parameters required for pCN: standard deviation and rho
#'@return list that contains updated states of the first and the second particles, updatd value of identical flag and the cost of calculations
#'@export

synchronous2_pcn_coupling <- function(chain_state1, chain_state2, identical, tuning){

  cost = 0  # cost of proposal generation
  # tuning parameters that define autoregressive proposal
  proposal_sd <- tuning$proposal_sd
  proposal_rho <- tuning$proposal_rho
  proposal_sd_factor <- sqrt(1-proposal_rho^2) * proposal_sd

  # sample common Gaussian increment
  noise <- rnorm(dimension)

  # sample proposals of fine chains
  state1 <- proposal_rho * chain_state1 + proposal_sd_factor * noise
  state2 <- proposal_rho * chain_state2 + proposal_sd_factor * noise

  # if chains are identical, they stay together;
  # if chains are not identical, the synchronous coupling cannot allow them to be

  return(list(state1 = state1, state2 = state2, identical = identical, cost = cost))

}

