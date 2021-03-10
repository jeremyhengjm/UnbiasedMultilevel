#' @rdname unbiased_increment
#' @title Unbiased estimator for the value of the increment at the given level
#' @description Generate four coupled chains and compute the estimator of the increment
#' @param level an integer that determines the target probability distributions of four chains
#' @param rinit function that is utilized for initialization of chains (for instance, samples from prior)
#' @param single_kernel a function that makes a single step through a specified MCMC kernel for a single chain
#' @param coupled2_kernel a function that makes a single step through a specified coupled MCMC kernel for 2 chains at different levels
#' @param coupled2_kernel a function that makes a single step through a specified coupled MCMC kernel for all four chains
#' @param proposal_coupling2 a function that generates proposal for a given input for two coupled chains (burn-in + lag)
#' @param proposal_coupling4 a function that generates proposal for a given input for four coupled chains (main run)
#' @param tuning a list of parameters reuired for MCMC iterations (for instance standard deviation for RWM)
#' @param tuning_coarse a list of parameters reuired for MCMC iterations for the coarse level (for instance standard deviation for RWM)
#' @param tuning_fine a list of parameters reuired for MCMC iterations for the fine level (for instance standard deviation for RWM)
#' @param h function that represents quantity of interest. Depends on level and spatial argument.
#' @param k an integer: lower bound for time-averaging
#' @param m an integer:  upper bound for time-averaging
#' @param max_iterations and integer: bound for the number of steps through coupled MCMC kernel
#' @param samping_factor a real value that controls the magnitude of perturnabtion at the initialization step
#'@return a list with the value of MCMC estimator without correction, value of Unbiased MCMC estimator, meeting time, value of iteration counter, flag that is "True" if chains have met before the iteration counter reached the value in max_iterations, cost of calculations
#'@export

unbiased_increment <- function(level, rinit, single_kernel, 
                               coupled2_kernel, coupled4_kernel, proposal_coupling2, proposal_coupling4,
                               tuning, tuning_coarse, tuning_fine,
                               h = function(l, x) x, k = 0, m = 1, max_iterations = Inf,
                               sampling_factor = 0.2){

  cost = 0  # number of operations
  # initialize chains
  state_coarse1 <- rinit(level-1)
  cost = cost + 2 ^ (level - 1)  # single calculation of the likelihood at level (l - 1)
  state_coarse2 <- rinit(level-1)
  cost = cost + 2 ^ (level - 1)  # single calculation of the likelihood at level (l - 1)
  state_fine1 <- rinit(level)
  cost = cost + 2 ^ level  # single calculation of the likelihood at level l
  state_fine2 <- rinit(level)
  cost = cost + 2 ^ level  # single calculation of the likelihood at level l
  identical_coarse <- FALSE
  identical_fine <- FALSE
  
  base_state = 0.0 + state_coarse1$chain_state
  generation <- TRUE
  while (generation)
  {
    pert = sampling_factor * rnorm(length(base_state))
    state_coarse1$chain_state <- base_state + pert
    state_coarse1$current_pdf <- logtarget(level, state_coarse1$chain_state)
    cost = cost + 2 ^ level
    if (is.finite(state_coarse1$current_pdf))
    {
      generation <- FALSE
    }
  }

  generation <- TRUE
  while (generation)
  {
    pert = sampling_factor * rnorm(length(base_state))
    state_fine1$chain_state <- base_state + pert
    state_fine1$current_pdf <- logtarget(level, state_fine1$chain_state)
    cost = cost + 2 ^ level
    if (is.finite(state_fine1$current_pdf))
    {
      generation <- FALSE
    }
  }

  base_state = 0.0 + state_coarse2$chain_state
  generation <- TRUE
  while (generation)
  {
    pert = sampling_factor * rnorm(length(base_state))
    state_coarse2$chain_state <- base_state + pert
    state_coarse2$current_pdf <- logtarget(level, state_coarse2$chain_state)
    cost = cost + 2 ^ level
    if (is.finite(state_coarse2$current_pdf))
    {
      generation <- FALSE
    }
  }

  generation <- TRUE
  while (generation)
  {
    pert = sampling_factor * rnorm(length(base_state))
    state_fine2$chain_state <- base_state + pert
    state_fine2$current_pdf <- logtarget(level, state_fine2$chain_state)
    cost = cost + 2 ^ level
    if (is.finite(state_fine2$current_pdf))
    {
      generation <- FALSE
    }
  }

  mcmcestimator_coarse <- h(level-1, state_coarse1$chain_state)
  mcmcestimator_fine <- h(level, state_fine1$chain_state)
  dimh <- length(mcmcestimator_coarse)
  if (k > 0){
    mcmcestimator_coarse <- rep(0, dimh)
    mcmcestimator_fine <- rep(0, dimh)
  }

  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction_coarse <- rep(0, dimh)
  correction_fine <- rep(0, dimh)
  # create the lag
  coupled_kernel_output <- coupled2_kernel(level,
                                           state_coarse1, state_fine1,
                                           identical = FALSE,
                                           tuning = tuning,
                                           proposal_coupling = proposal_coupling2)
  state_coarse1 <- coupled_kernel_output$state1
  state_fine1 <- coupled_kernel_output$state2
  cost = cost + coupled_kernel_output$cost

  # coupled_kernel_output <- coupled2_kernel(level,
  #                                          state_coarse1, state_fine1,
  #                                          identical = FALSE,
  #                                          tuning = tuning_coarse,
  #                                          proposal_coupling = proposal_coupling2)
  if (k == 0){
    correction_coarse <- correction_coarse + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h(level-1, state_coarse1$chain_state) - h(level-1, state_coarse2$chain_state))
    correction_fine <- correction_fine + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h(level, state_fine1$chain_state) - h(level, state_fine2$chain_state))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator_coarse <- mcmcestimator_coarse + h(level-1, state_coarse1$chain_state)
    mcmcestimator_fine <- mcmcestimator_fine + h(level, state_fine1$chain_state)
  }

  # initialize
  iter <- 1
  meet_coarse <- FALSE
  meet_fine <- FALSE
  finished <- FALSE
  meetingtime_coarse <- Inf
  meetingtime_fine <- Inf

  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    # increment counter
    iter <- iter + 1

    # run coupled kernel
    coupled_kernel_output <- coupled4_kernel(level,
                                             state_coarse1, state_coarse2,
                                             state_fine1, state_fine2,
                                             identical_coarse, identical_fine,
                                             tuning, tuning,
                                             proposal_coupling4)
  
    state_coarse1 <- coupled_kernel_output$state_coarse1
    state_coarse2 <- coupled_kernel_output$state_coarse2
    state_fine1 <- coupled_kernel_output$state_fine1
    state_fine2 <- coupled_kernel_output$state_fine2
    identical_coarse <- coupled_kernel_output$identical_coarse
    identical_fine <- coupled_kernel_output$identical_fine
    cost = cost + coupled_kernel_output$cost

    # update estimator for coarse discretization level
    if (meet_coarse){
      if (k <= iter && iter <= m){
        mcmcestimator_coarse <- mcmcestimator_coarse + h(level-1, state_coarse1$chain_state)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_coarse <- mcmcestimator_coarse + h(level-1, state_coarse1$chain_state)
        }
        correction_coarse <- correction_coarse + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h(level-1, state_coarse1$chain_state) - h(level-1, state_coarse2$chain_state))
      }
    }

    # update estimator for fine discretization level
    if (meet_fine){
      if (k <= iter && iter <= m){
        mcmcestimator_fine <- mcmcestimator_fine + h(level, state_fine1$chain_state)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_fine <- mcmcestimator_fine + h(level, state_fine1$chain_state)
        }
        correction_fine <- correction_fine + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h(level, state_fine1$chain_state) - h(level, state_fine2$chain_state))
      }
    }

    # check if meeting occurs for coarse discretization level
    if (identical_coarse && !meet_coarse){
      meet_coarse <- TRUE # recording meeting time tau_coarse
      meetingtime_coarse <- iter
    }

    # check if meeting occurs for fine discretization level
    if (identical_fine && !meet_fine){
      meet_fine <- TRUE # recording meeting time tau_fine
      meetingtime_fine <- iter
    }

    # stop after max(m, tau_coarse, tau_fine) steps
    if (iter >= max(meetingtime_coarse, meetingtime_fine, m)){
      finished <- TRUE
    }
  }
  # compute mcmc estimators and their difference
  mcmcestimator_coarse <- mcmcestimator_coarse / (m - k + 1)
  mcmcestimator_fine <- mcmcestimator_fine / (m - k + 1)
  mcmcestimator <- mcmcestimator_fine - mcmcestimator_coarse

  # compute unbiased estimators and their difference
  uestimator_coarse <- mcmcestimator_coarse + correction_coarse
  uestimator_fine <- mcmcestimator_fine + correction_fine
  uestimator <- uestimator_fine - uestimator_coarse

  return(list(mcmcestimator_coarse = mcmcestimator_coarse, mcmcestimator_fine = mcmcestimator_fine,
              uestimator_coarse = uestimator_coarse, uestimator_fine = uestimator_fine,
              mcmcestimator = mcmcestimator, uestimator = uestimator,
              meetingtime_coarse = meetingtime_coarse, meetingtime_fine = meetingtime_fine,
              iteration = iter, finished = finished, cost = cost))

}

