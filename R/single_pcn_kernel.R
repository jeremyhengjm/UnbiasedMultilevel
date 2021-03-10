#' @rdname single_pcn_kernel
#' @title pCN kernel for a single MCMC chain that target probability distribution at the given level
#' @description Generation of pCN proposal and accept/reject step
#' @param level a integer that determines density of probability distribution inthe multi-level approach
#' @param state a list with current position of the particle (element of the chain) and log of the density of probability distribution
#' @param tuning a list of parameters for pCN iteration: standard deviation and rho
#'@return a updated state of the chain and cost of computations
#'@export

single_pcn_kernel <- function(level, state, tuning){
  # extract state and pdf
  chain_state <- state$chain_state
  current_pdf <- state$current_pdf

  # tuning parameters that define autoregressive proposal
  proposal_sd <- tuning$proposal_sd
  proposal_rho <- tuning$proposal_rho
  proposal_sd_factor <- sqrt(1-proposal_rho^2) * proposal_sd

  # sample proposal and compute pdf
  proposal_value <- proposal_rho * chain_state + proposal_sd_factor * rnorm(dimension)
  proposal_pdf <- logtarget(level, proposal_value)

  if(is.finite(proposal_pdf))
  {
    # compute acceptance probability
     logacceptprob <- proposal_pdf - current_pdf +
       sum(dnorm(chain_state, mean = proposal_rho * proposal_value, sd = proposal_sd_factor, log = TRUE)) -
       sum(dnorm(proposal_value, mean = proposal_rho * chain_state, sd = proposal_sd_factor, log = TRUE))

     # accept or reject proposal
     accept <- (log(runif(1)) < logacceptprob)
  }
  else
  {
     accept <- FALSE
  }
  if (accept){
    chain_data = list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept)    
  } else {
    chain_data = list(chain_state = chain_state, current_pdf = current_pdf, accept = accept)
  }
  return(list(chain_data = chain_data, cost = 2 ^ level))
}

