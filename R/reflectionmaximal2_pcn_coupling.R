#' @rdname reflectionmaximal2_pcn_coupling
#' @title Reflection maximal coupling of 2 pCN proposals
#' @description Generate pCN proposals for two chains via reflection maximal coupling
#' @param chain_state1 a vector with coordinates of the first particle
#' @param chain_state2 a vector with coordinates of the second particle
#' @param identical a flag that is True if chains are identical and False otherwise
#' @param tuning a list that contains parameters needed for pCN: standard devation and rho
#'@return a list that contains state of the first chain, state of the second chain, updated value of the flag "identical", cost of proposal generation
#'@export

reflectionmaximal2_pcn_coupling <- function(chain_state1, chain_state2, identical, tuning){
  cost <- 0  # cost of proposal generation
  # extract tuning parameters
  proposal_sd <- tuning$proposal_sd
  proposal_rho <- tuning$proposal_rho
  proposal_sd_factor <- sqrt(1-proposal_rho^2) * proposal_sd

  # sample first proposal
  randn1 <- rnorm(dimension)
  state1 <- proposal_rho * chain_state1 + proposal_sd_factor * randn1

  # difference
  zdiff <- proposal_rho * (chain_state1 - chain_state2) / proposal_sd_factor

  # evaluate proposal transition densities at first proposal
  pdf1 <- sum(dnorm(randn1, log = TRUE))
  pdf2 <- sum(dnorm(randn1 + zdiff, log = TRUE))
  logacceptprob <- min(pdf1, pdf2) - pdf1

  if (log(runif(1)) < logacceptprob){
    # return common proposal for both chains
    return(list(state1 = state1, state2 = state1, identical = TRUE, cost = cost))

  } else {
    # perform reflection for second proposal
    evec <- zdiff / sqrt(sum(zdiff^2))
    randn2 <- randn1 - 2 * sum(evec * randn1) * evec
    state2 <- proposal_rho * chain_state2 + proposal_sd_factor * randn2
    return(list(state1 = state1, state2 = state2, identical = FALSE, cost = cost))

  }
}

