#' @rdname synchronous2_hmc_coupling_extra_level
#' @title Syncronous coupling of two HMC proposals for chains that target distributions at different levels
#' @description Sampling of common momentum for both particles and leapfrog integration at each level
#' @param level an iteger number that determines probability distribution in the multi-level approach
#' @param chain_state1 a vector with the position of the first particle (at the "level" (level - 1))
#' @param chain_state2 a vector with the position of the second particle (at the "level" level)
#' @param identical a boolean variable takes the value True if two chains coincide and False otherwise. This parameter is always False and is left for compatibility with other functions in the package
#' @param tuning a list that contains parameters required for leapfrog integration: stepsize and number os steps
#'@return a list with updated values of the states of the first and the second chain, velocites of each of the particles, initial velocity (random sample from normal distribution), updated value of the "identical" flag and cost of computations
#'@export

synchronous2_hmc_coupling_extra_level <- function(level, chain_state1, chain_state2, identical, tuning){

  cost = 0  # cost of proposal generation
  # sample common velocity or momentum
  initial_velocity <- rnorm(dimension)

  # run leapfrog integrator
  leapfrog_result <- leapfrog(level - 1, chain_state1, initial_velocity, tuning)
  state1 <- leapfrog_result$x
  velocity1 <- - leapfrog_result$v
  cost = cost + leapfrog_result$cost

  leapfrog_result <- leapfrog(level, chain_state2, initial_velocity, tuning)
  state2 <- leapfrog_result$x
  velocity2 <- - leapfrog_result$v
  cost = cost + leapfrog_result$cost

  # if chains are identical, they stay together;
  # if chains are not identical, the synchronous coupling cannot allow them to be

  return(list(initial_velocity1 = initial_velocity, initial_velocity2 = initial_velocity,
              state1 = state1, state2 = state2, velocity1 = velocity1, velocity2 = velocity2,
              identical = identical, cost = cost))
}

