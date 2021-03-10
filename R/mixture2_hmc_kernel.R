#' @rdname mixture2_hmc_kernel
#' @title Mixture of kernels with HMC synchronous coupling and RWM maximal coupling for two coupled chains at the given level
#' @description Either HMC kernel or RWM kernel is selected for MCMC iteration with a certain probability
#' @param level an integer that controls the probability distribution in a multi-level setting
#' @param state1 a list that contains vector of coordinates and log of the density of probability distribution for the first chain
#' @param state2 a list that contains vector of coordinates and log of the density of probability distribution for the second chain
#' @param identical a boolean variable: "True" if chains coincide, "False" otherwise
#' @param tuning list of parameters for HMC and RWM kernels (number os steps and stepsize for HMC, standard devation for RWM)
#'@return list that contains updated states of the first and the second chain, updated value of the flag "identical" and cost of computations
#'@export

mixture2_hmc_kernel <- function(level, state1, state2, identical, tuning, proposal_coupling){
  if (runif(1) < probability_maximal_coupling){
    return(coupled2_rwmh_kernel(level, state1, state2, identical, tuning, proposal_coupling))
  } else {
    return(coupled2_hmc_kernel(level, state1, state2, identical, tuning, synchronous2_hmc_coupling))
  }
}

