#' @rdname unbiased_expectation
#' @title Unbiased estimator for the value at level zero
#' @description Generate two coupled chains and compute the estimator
#' @param level an integer that determines the target probability distribution of two chains
#' @param rinit function that is utilized for initialization of chains (for instance, samples from prior)
#' @param single_kernel function that makes a single step through a specified MCMC kernel
#' @param tuning a list of parameters reuired for MCMC iterations (for instance standard deviation for RWM)
#' @param coupled_kernel function that makes a single step through a specified coupled MCMC kernel (for 2 chains)
#' @param proposal_coupling function that generates proposal for a given input
#' @param h function that represents quantity of interest. Depends on level and spatial argument.
#' @param k an integer: lower bound for time-averaging
#' @param m an integer:  upper bound for time-averaging
#' @param max_iterations and integer: bound for the number of steps through coupled MCMC kernel
#'@return a list with the value of MCMC estimator without correction, value of Unbiased MCMC estimator, meeting time, value of iteration counter, flag that is "True" if chains have met before the iteration counter reached the value in max_iterations, cost of calculations
#'@export

unbiased_expectation <- function(level, rinit, single_kernel, tuning,
                                 coupled_kernel, proposal_coupling,
                                 h = function(l, x) x, k = 0, m = 1, max_iterations = Inf){

  cost = 0  # number of operations
  # initialize chains
  state1 <- rinit(level)
  cost = cost + 2 ^ level  # single evaluation of loglikelihood
  state2 <- rinit(level)
  cost = cost + 2 ^ level  # single evaluation of loglikelihood
  identical = FALSE

  #coupled_kernel_output = coupled_kernel(level, state1, state2, identical, tuning, proposal_coupling)
  #identical <- FALSE
  #state2 <- coupled_kernel_output$state1

  # initialize estimator computation
  mcmcestimator <- h(level, state1$chain_state)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }

  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  single_kernel_output = single_kernel(level, state1, tuning)
  state1 <- single_kernel_output$chain_data
  cost = cost + single_kernel_output$cost
  # I have to go to tuning directly and return number of function evaluation
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h(level, state1$chain_state) - h(level, state2$chain_state))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(level, state1$chain_state)
  }

  # initialize
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf

  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    # increment counter
    iter <- iter + 1

    # run coupled kernel
    coupled_kernel_output <- coupled_kernel(level,
                                            state1, state2,
                                            identical,
                                            tuning,
                                            proposal_coupling)

    state1 <- coupled_kernel_output$state1
    state2 <- coupled_kernel_output$state2
    identical <- coupled_kernel_output$identical
    cost = cost + coupled_kernel_output$cost
    #print("chain_state1")
    #print(state1)
    #print("chain_state2")
    #print(state2)
    #print("chain_diff")
    #print(state2$chain_state - state1$chain_state)
    #print("chain_state2")
    #print(state2)

    # update estimator for fine discretization level
    if (meet){
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(level, state1$chain_state)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(level, state1$chain_state)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h(level, state1$chain_state) - h(level, state2$chain_state))
      }
    }

    # check if meeting occurs
    if (identical && !meet){
      meet <- TRUE # recording meeting time
      meetingtime <- iter
    }

    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  # compute mcmc estimators
  mcmcestimator <- mcmcestimator / (m - k + 1)
  #print("mcmcestimator")
  #print(mcmcestimator)
  #print("correction")
  #sprint(correction)

  # compute unbiased estimator
  uestimator <- mcmcestimator + correction

  return(list(mcmcestimator = mcmcestimator, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished, cost = cost))

}

