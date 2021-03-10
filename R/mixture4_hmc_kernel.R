#' @rdname mixture4_hmc_kernel
#' @title Mixture of kernels with HMC synchronous coupling and RWM maximal coupling for four coupled chains
#' @description Either HMC kernel or RWM kernel is selected for MCMC iteration with a certain probability
#' @param level an integer that controls the probability distribution in a multi-level setting
#' @param state_coarse1 a list that contains vector of coordinates and log of the density of probability distribution for the first chain at the coarse level: "level - 1"
#' @param state_coarse2 a list that contains vector of coordinates and log of the density of probability distribution for the second chain at the coarse level: "level - 1"
#' @param state_fine1 a list that contains vector of coordinates and log of the density of probability distribution for the first chain at the fine level: "level"
#' @param state_fine2 a list that contains vector of coordinates and log of the density of probability distribution for the second chain at the fine level: "level"
#' @param identical_coarse a boolean variable: "True" if chains at the coarse level coincide, "False" otherwise
#' @param identical_fine a boolean variable: "True" if chains at the fine level coincide, "False" otherwise
#' @param tuning_coarse list of parameters for HMC and RWM kernels (number os steps and stepsize for HMC, standard devation for RWM) for two chains at the coarse level
#' @param tuning_fine list of parameters for HMC and RWM kernels (number os steps and stepsize for HMC, standard devation for RWM) for two chains at the fine level
#'@return list that contains updated states of the first and the second chains at the coarse level,  updated states of the first and the second chains at the fine level, updated value of the flag "identical_coarse", updated value of the flag "identical_fine" and cost of computations
#'@export

mixture4_hmc_kernel <- function(level,
                                state_coarse1, state_coarse2,
                                state_fine1, state_fine2,
                                identical_coarse, identical_fine,
                                tuning_coarse, tuning_fine,
                                proposal_coupling){

  if (runif(1) < probability_maximal_coupling){
    return(coupled4_rwmh_kernel(level,
                                state_coarse1, state_coarse2,
                                state_fine1, state_fine2,
                                identical_coarse, identical_fine,
                                tuning_coarse, tuning_fine,
                                proposal_coupling))
  } else {
    return(coupled4_hmc_kernel(level,
                               state_coarse1, state_coarse2,
                               state_fine1, state_fine2,
                               identical_coarse, identical_fine,
                               tuning_coarse, tuning_fine,
                               synchronous4_hmc_coupling))
  }

}
