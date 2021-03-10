#' @rdname single_mixture_hmc_kernel
#' @title Mixture of kernels with HMC synchronous coupling and RWM maximal coupling a single chain at the given level
#' @description Either HMC kernel or RWM kernel is selected for MCMC iteration with a certain probability
#' @param level an integer that controls the probability distribution in a multi-level setting
#' @param state a list that contains vector of coordinates and log of the density of probability distribution for the first chain
#' @param tuning list of parameters for HMC and RWM kernels (number os steps and stepsize for HMC, standard devation for RWM)
#'@return list that contains updated states of the first and the second chain and cost of computations
#'@export

single_mixture_hmc_kernel <- function(level, state, tuning){
  if (runif(1) < probability_maximal_coupling){
    return(single_rwmh_kernel(level, state, tuning))
  } else{
    return(single_hmc_kernel(level, state, tuning))
  }
}

