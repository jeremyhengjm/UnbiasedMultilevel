#' @rdname coxprocess_unbiased_expectation
#' @title Unbiased estimator for the value at level zero
#' @description Generate two coupled chains and compute the estimator
#' @param parameters a list of model parameters (sigmasq, mu, beta)
#' @param model a list of model functions with keys \code{logtarget}, \code{gradlogtarget}, \code{rinit} and \code{dimension}
#' @param mcmc a list of MCMC kernels with keys \code{single_kernel} and \code{coupled2_kernel}
#' @param h function that represents quantity of interest
#' @param k an integer: lower bound for time-averaging
#' @param m an integer:  upper bound for time-averaging
#' @param max_iterations iteration at which to stop the while loop (default to infinity)
#'@return a list with the value of MCMC estimator without correction, value of Unbiased MCMC estimator, meeting time, value of iteration counter, flag that is "True" if chains have met before the iteration counter reached the value in max_iterations, cost of calculations
#'@export

coxprocess_unbiased_expectation <- function(parameters, model, mcmc, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  # MCMC kernels
  single_kernel <- mcmc$single_kernel
  coupled2_kernel <- mcmc$coupled2_kernel

  # initialize chains
  chain_state1 <- model$rinit(rnorm(model$dimension))
  chain_state2 <- model$rinit(rnorm(model$dimension))
  current_pdf1 <- model$logtarget(chain_state1)
  current_pdf2 <- model$logtarget(chain_state2)
  identical <- FALSE

  # initialize estimator computation
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }

  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  single_kernel_output <- single_kernel(chain_state1, current_pdf1)
  chain_state1 <- single_kernel_output$chain_state
  current_pdf1 <- single_kernel_output$current_pdf
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
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
    coupled2_kernel_output <- coupled2_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
    chain_state1 <- coupled2_kernel_output$chain_state1
    chain_state2 <- coupled2_kernel_output$chain_state2
    current_pdf1 <- coupled2_kernel_output$current_pdf1
    current_pdf2 <- coupled2_kernel_output$current_pdf2
    identical <- all(chain_state1 == chain_state2)

    # update estimator for fine discretization level
    if (meet){
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h(chain_state1) - h(chain_state2))
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

  # compute unbiased estimator
  unbiasedestimator <- mcmcestimator + correction

  # compute cost in units of single kernel at current discretization level
  cost <- 2 * (meetingtime - 1) + max(1, m + 1 - meetingtime)

  return(list(mcmcestimator = mcmcestimator, unbiasedestimator = unbiasedestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished, cost = cost))

}

#' @rdname coxprocess_unbiased_increment
#' @title Unbiased estimator for the value of the increment at the given level
#' @description Generate four coupled chains and compute the estimator of the increment
#' @param parameters a list of model parameters (sigmasq, mu, beta)
#' @param model_coarse a list of functions for coarser resolution model with keys \code{logtarget}, \code{gradlogtarget}, \code{rinit} and \code{dimension}
#' @param model_fine a list of functions for finer resolution model with keys \code{logtarget}, \code{gradlogtarget}, \code{rinit} and \code{dimension}
#' @param mcmc a list of MCMC kernels with keys \code{multilevel_kernel} and \code{coupled4_kernel}
#' @param h_coarse function that represents quantity of interest on coarser resolution
#' @param h_fine function that represents quantity of interest on finer resolution
#' @param k an integer: lower bound for time-averaging
#' @param m an integer:  upper bound for time-averaging
#' @param max_iterations iteration at which to stop the while loop (default to infinity)
#'@return a list with the value of MCMC estimator without correction, value of Unbiased MCMC estimator, meeting time, value of iteration counter, flag that is "True" if chains have met before the iteration counter reached the value in max_iterations, cost of calculations
#'@export

coxprocess_unbiased_increment <- function(parameters, model_coarse, model_fine, mcmc,
                                          h_coarse = function(x) x, h_fine = function(x) x,
                                          k = 0, m = 1, max_iterations = Inf){
  # MCMC kernels
  multilevel_kernel <- mcmc$multilevel_kernel
  coupled4_kernel <- mcmc$coupled4_kernel
  rnorm_coarse <- mcmc$rnorm_coarse

  # initialize chains
  rnorm_fine1 <- rnorm(model_fine$dimension)
  rnorm_fine2 <- rnorm(model_fine$dimension)
  chain_state_fine1 <- model_fine$rinit(rnorm_fine1)
  chain_state_fine2 <- model_fine$rinit(rnorm_fine2)

  rnorm_coarse1 <- rnorm_coarse(rnorm_fine1)
  rnorm_coarse2 <- rnorm_coarse(rnorm_fine2)
  chain_state_coarse1 <- model_coarse$rinit(rnorm_coarse1)
  chain_state_coarse2 <- model_coarse$rinit(rnorm_coarse2)

  current_pdf_coarse1 <- model_coarse$logtarget(chain_state_coarse1)
  current_pdf_coarse2 <- model_coarse$logtarget(chain_state_coarse2)
  current_pdf_fine1 <- model_fine$logtarget(chain_state_fine1)
  current_pdf_fine2 <- model_fine$logtarget(chain_state_fine2)

  identical_coarse <- FALSE
  identical_fine <- FALSE

  mcmcestimator_coarse <- h_coarse(chain_state_coarse1)
  mcmcestimator_fine <- h_fine(chain_state_fine1)
  dimh <- length(mcmcestimator_coarse)
  if (k > 0){
    mcmcestimator_coarse <- rep(0, dimh)
    mcmcestimator_fine <- rep(0, dimh)
  }

  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction_coarse <- rep(0, dimh)
  correction_fine <- rep(0, dimh)

  multilevel_kernel_output <- multilevel_kernel(chain_state_coarse1, current_pdf_coarse1, chain_state_fine1, current_pdf_fine1)
  chain_state_coarse1 <- multilevel_kernel_output$chain_state_coarse
  current_pdf_coarse1 <- multilevel_kernel_output$current_pdf_coarse
  chain_state_fine1 <- multilevel_kernel_output$chain_state_fine
  current_pdf_fine1 <- multilevel_kernel_output$current_pdf_fine

  if (k == 0){
    correction_coarse <- correction_coarse + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h_coarse(chain_state_coarse1) - h_coarse(chain_state_coarse2))
    correction_fine <- correction_fine + (min(1, (0 - k + 1)/(m - k + 1))) *
      (h_fine(chain_state_fine1) - h_fine(chain_state_fine2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator_coarse <- mcmcestimator_coarse + h_coarse(chain_state_coarse1)
    mcmcestimator_fine <- mcmcestimator_fine + h_fine(chain_state_fine1)
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
    coupled4_kernel_output <- coupled4_kernel(chain_state_coarse1, chain_state_coarse2,
                                              current_pdf_coarse1, current_pdf_coarse2,
                                              chain_state_fine1, chain_state_fine2,
                                              current_pdf_fine1, current_pdf_fine2)
    chain_state_coarse1 <- coupled4_kernel_output$chain_state_coarse1
    chain_state_coarse2 <- coupled4_kernel_output$chain_state_coarse2
    current_pdf_coarse1 <- coupled4_kernel_output$current_pdf_coarse1
    current_pdf_coarse2 <- coupled4_kernel_output$current_pdf_coarse2

    chain_state_fine1 <- coupled4_kernel_output$chain_state_fine1
    chain_state_fine2 <- coupled4_kernel_output$chain_state_fine2
    current_pdf_fine1 <- coupled4_kernel_output$current_pdf_fine1
    current_pdf_fine2 <- coupled4_kernel_output$current_pdf_fine2

    identical_coarse <- all(chain_state_coarse1 == chain_state_coarse2)
    identical_fine <- all(chain_state_fine1 == chain_state_fine2)

    # update estimator for coarse discretization level
    if (meet_coarse){
      if (k <= iter && iter <= m){
        mcmcestimator_coarse <- mcmcestimator_coarse + h_coarse(chain_state_coarse1)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_coarse <- mcmcestimator_coarse + h_coarse(chain_state_coarse1)
        }
        correction_coarse <- correction_coarse + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h_coarse(chain_state_coarse1) - h_coarse(chain_state_coarse2))
      }
    }

    # update estimator for fine discretization level
    if (meet_fine){
      if (k <= iter && iter <= m){
        mcmcestimator_fine <- mcmcestimator_fine + h_fine(chain_state_fine1)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_fine <- mcmcestimator_fine + h_fine(chain_state_fine1)
        }
        correction_fine <- correction_fine + (min(1, (iter-1 - k + 1)/(m - k + 1))) *
          (h_fine(chain_state_fine1) - h_fine(chain_state_fine2))
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
  unbiasedestimator_coarse <- mcmcestimator_coarse + correction_coarse
  unbiasedestimator_fine <- mcmcestimator_fine + correction_fine
  unbiasedestimator <- unbiasedestimator_fine - unbiasedestimator_coarse

  # compute cost in units of CPF kernel at current discretization level
  cost_coarse <- 2 * (meetingtime_coarse - 1) + max(1, m + 1 - meetingtime_coarse)
  cost_fine <- 2 * (meetingtime_fine - 1) + max(1, m + 1 - meetingtime_fine)


  return(list(mcmcestimator_coarse = mcmcestimator_coarse, mcmcestimator_fine = mcmcestimator_fine,
              unbiasedestimator_coarse = unbiasedestimator_coarse, unbiasedestimator_fine = unbiasedestimator_fine,
              mcmcestimator = mcmcestimator, unbiasedestimator = unbiasedestimator,
              meetingtime_coarse = meetingtime_coarse, meetingtime_fine = meetingtime_fine,
              iteration = iter, finished = finished,
              cost_coarse = cost_coarse, cost_fine = cost_fine))

}






