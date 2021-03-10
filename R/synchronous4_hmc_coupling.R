#' @rdname synchronous4_hmc_coupling
#' @title Syncronous coupling of four HMC proposals
#' @description Sampling of common momentum for all four particles and leapfrog integration
#' @param level an integer number that determines probability distribution in the multi-level approach
#' @param chain_state_coarse1 a vector with the position of the first particle at the coarse level (level - 1)
#' @param chain_state_coarse2 a vector with the position of the second particle at the coarse level (level - 1)
#' @param chain_state_fine1 a vector with the position of the first particle at the fine level (level)
#' @param chain_state_fine2 a vector with the position of the second particle at the fine level (level)
#' @param identical_coarse a boolean variable takes the value True if two chains at the coarse level coincide and False otherwise
#' @param identical_fine a boolean variable takes the value True if two chains at the fine level coincide and False otherwise
#' @param tuning_coarse a list that contains parameters required for leapfrog integration at the coarse level: stepsize and number os steps
#' @param tuning_fine a list that contains parameters required for leapfrog integration at the fine level: stepsize and number os steps
#'@return a list with updated values of states of the first and the second chain at the coarse level, updated values of states of the first and the second chain at the fine level, velocites of each of the particles at the coarse level, velocites of each of the particles at the fine level, initial velocity (random sample from normal distribution), updated value of the "identical_coarse" flag, updated value of the "identical_fine" flag and cost of computations
#'@export

synchronous4_hmc_coupling <- function(level,
                                      chain_state_coarse1, chain_state_coarse2,
                                      chain_state_fine1, chain_state_fine2,
                                      identical_coarse, identical_fine,
                                      tuning_coarse, tuning_fine){

  cost = 0  # cost of proposal generation
  # sample common velocity or momentum
  initial_velocity <- rnorm(dimension)

  # run leapfrog integrator for coarse level
  leapfrog_result <- leapfrog(level-1, chain_state_coarse1, initial_velocity, tuning_coarse)
  state_coarse1 <- leapfrog_result$x
  velocity_coarse1 <- - leapfrog_result$v
  cost = cost + leapfrog_result$cost
  if (identical_coarse){
    state_coarse2 <- state_coarse1
    velocity_coarse2 <- velocity_coarse1
  } else {
    leapfrog_result <- leapfrog(level-1, chain_state_coarse2, initial_velocity, tuning_coarse)
    state_coarse2 <- leapfrog_result$x
    velocity_coarse2 <- - leapfrog_result$v
    cost = cost + leapfrog_result$cost
  }

  # run leapfrog integrator for fine level
  leapfrog_result <- leapfrog(level, chain_state_fine1, initial_velocity, tuning_fine)
  state_fine1 <- leapfrog_result$x
  velocity_fine1 <- - leapfrog_result$v
  cost = cost + leapfrog_result$cost
  if (identical_fine){
    state_fine2 <- state_fine1
    velocity_fine2 <- velocity_fine1
  } else {
    leapfrog_result <- leapfrog(level, chain_state_fine2, initial_velocity, tuning_fine)
    state_fine2 <- leapfrog_result$x
    velocity_fine2 <- - leapfrog_result$v
    cost = cost + leapfrog_result$cost
  }

  # if chains are identical, they stay together;
  # if chains are not identical, the synchronous coupling cannot allow them to be

  return(list(initial_velocity_coarse1 = initial_velocity, initial_velocity_coarse2 = initial_velocity,
              state_coarse1 = state_coarse1, state_coarse2 = state_coarse2,
              velocity_coarse1 = velocity_coarse1, velocity_coarse2 = velocity_coarse2,
              initial_velocity_fine1 = initial_velocity, initial_velocity_fine2 = initial_velocity,
              state_fine1 = state_fine1, state_fine2 = state_fine2,
              velocity_fine1 = velocity_fine1, velocity_fine2 = velocity_fine2,
              identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))
}
