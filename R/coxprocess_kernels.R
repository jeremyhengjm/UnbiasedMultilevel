#'@rdname coxprocess_kernels
#'@title Markov chain Monte Carlo kernels for log-Gaussian Cox process model
#'@description This function defines Markov chain Monte Carlo kernels and its couplings
#'@param model list containing \code{logtarget}, \code{gradlogtarget} and \code{dimension}
#'@param tuning list containing \code{stepsize} and \code{nsteps}
#'@return a list containing the keys \code{single_hmc_kernel}, \code{coupled2_hmc_kernel}, \code{single_kernel} and \code{coupled2_kernel}
#'@export
coxprocess_kernels <- function(model, tuning){

  # get log-posterior density and its gradient
  logtarget <- model$logtarget
  gradlogtarget <- model$gradlogtarget
  dimension <- model$dimension

  # tuning parameters
  stepsize <- tuning$stepsize # step size in the leap-frog integrator
  nsteps <- tuning$nsteps # number of leap-frog steps
  proposal_sd <- tuning$proposal_sd # proposal standard deviation of RWMH kernel
  probability_rwmh <- tuning$probability_rwmh # probability of selecting RWMH kernel in a mixture

  # leap frog integrator
  leapfrog_integrator <- function(x, v){
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  # single HMC kernel
  single_hmc_kernel <- function(chain_state, current_pdf){
    current_v <- rnorm(dimension) # velocity or momentum
    leapfrog_result <- leapfrog_integrator(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x
    proposed_pdf <- logtarget(proposed_x)
    accept_ratio <- proposed_pdf - current_pdf
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
    accept <- FALSE
    if (is.finite(accept_ratio)){
      accept <- (log(runif(1)) < accept_ratio)
    }

    if (accept){
      chain_state <- proposed_x
      current_pdf <- proposed_pdf
      accept <- TRUE
    }
    return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
  }

  # coupled2 HMC kernel
  coupled2_hmc_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    current_v <- rnorm(dimension) # velocity or momentum, shared by both chains
    leapfrog_result1 <- leapfrog_integrator(chain_state1, current_v)
    leapfrog_result2 <- leapfrog_integrator(chain_state2, current_v)
    proposed_v1 <- - leapfrog_result1$v
    proposed_x1 <- leapfrog_result1$x
    proposed_v2 <- - leapfrog_result2$v
    proposed_x2 <- leapfrog_result2$x

    proposed_pdf1 <- logtarget(proposed_x1)
    proposed_pdf2 <- logtarget(proposed_x2)
    accept_ratio1 <- proposed_pdf1 - current_pdf1
    accept_ratio2 <- proposed_pdf2 - current_pdf2
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    current_kinetic_energy <- sum(current_v^2) / 2
    accept_ratio1 <- accept_ratio1 + current_kinetic_energy - sum(proposed_v1^2) / 2
    accept_ratio2 <- accept_ratio2 + current_kinetic_energy - sum(proposed_v2^2) / 2
    logu <- log(runif(1)) # shared by both chains
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(accept_ratio1)){
      accept1 <- (logu < accept_ratio1)
    }

    if (accept1){
      chain_state1 <- proposed_x1
      current_pdf1 <- proposed_pdf1
    }

    if (is.finite(accept_ratio2)){
      accept2 <- (logu < accept_ratio2)
    }

    if (accept2){
      chain_state2 <- proposed_x2
      current_pdf2 <- proposed_pdf2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
                accept1 = accept1, accept2 = accept2))
  }

  # single RWMH kernel
  single_rwmh_kernel <- function(chain_state, current_pdf){
    proposal_value <- chain_state + proposal_sd * rnorm(dimension)
    proposal_pdf <- logtarget(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
    }
  }

  reflectionmaximal2_coupling <- function(chain_state1, chain_state2){
    # sample first proposal
    randn1 <- rnorm(dimension)
    state1 <- chain_state1 + proposal_sd * randn1

    # difference
    zdiff <- (chain_state1 - chain_state2) / proposal_sd

    # evaluate proposal transition densities at first proposal
    pdf1 <- sum(dnorm(randn1, log = TRUE))
    pdf2 <- sum(dnorm(randn1 + zdiff, log = TRUE))
    logacceptprob <- min(pdf1, pdf2) - pdf1

    if (log(runif(1)) < logacceptprob){
      # return common proposal for both chains
      return(list(state1 = state1, state2 = state1))

    } else {
      # perform reflection for second proposal
      evec <- zdiff / sqrt(sum(zdiff^2))
      randn2 <- randn1 - 2 * sum(evec * randn1) * evec
      state2 <- chain_state2 + proposal_sd * randn2

      return(list(state1 = state1, state2 = state2))
    }
  }


  # coupled2 RWMH kernel
  coupled2_rwmh_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    # sample from maximal coupling
    proposal_value <- reflectionmaximal2_coupling(chain_state1, chain_state2)
    proposal1 <- proposal_value$state1
    proposal2 <- proposal_value$state2
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)

    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
                accept1 = accept1, accept2 = accept2))
  }

  # mixture of HMC and RWMH kernels
  single_mixture_kernel <- function(chain_state, current_pdf){
    if (runif(1) < probability_rwmh){
      return(single_rwmh_kernel(chain_state, current_pdf))
    } else{
      return(single_hmc_kernel(chain_state, current_pdf))
    }
  }

  # mixture of coupled2 HMC and RWMH kernels
  mixture_coupled2_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    if (runif(1) < probability_rwmh){
      return(coupled2_rwmh_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2))
    } else {
      return(coupled2_hmc_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2))
    }
  }

  return(list(single_hmc_kernel = single_hmc_kernel, coupled2_hmc_kernel = coupled2_hmc_kernel,
              single_kernel = single_mixture_kernel, coupled2_kernel = mixture_coupled2_kernel))
}

#'@rdname coxprocess_multilevel_kernels
#'@title Multilevel Markov chain Monte Carlo kernels for log-Gaussian Cox process model
#'@description This function defines multilevel Markov chain Monte Carlo kernels
#'@param model_coarse list for coarser resolution model containing \code{logtarget}, \code{gradlogtarget} and \code{dimension}
#'@param model_finer list for finer resolution model containing \code{logtarget}, \code{gradlogtarget} and \code{dimension}
#'@param tuning list containing \code{stepsize} and \code{nsteps}
#'@return a list containing the keys \code{single_kernel}, \code{coupled2_kernel}
#'@export
coxprocess_multilevel_kernels <- function(model_coarse, model_fine, tuning){

  # tuning parameters
  stepsize <- tuning$stepsize # step size in the leap-frog integrator
  nsteps <- tuning$nsteps # number of leap-frog steps
  proposal_sd <- tuning$proposal_sd # proposal standard deviation of RWMH kernel
  probability_rwmh <- tuning$probability_rwmh # probability of selecting RWMH kernel in a mixture

  # leap frog integrator
  leapfrog_integrator <- function(model, x, v){
    v <- v + stepsize * model$gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      if (step != nsteps){
        v <- v + stepsize * model$gradlogtarget(x)
      }
    }
    v <- v + stepsize * model$gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  # obtain standard Normal sample for coarser resolution by combining sample for finer resolution
  index_start <- seq(from = 1, to = model_fine$dimension-1, by = 2 * model_fine$ngrid)
  index_end <- seq(from = 2*model_fine$ngrid, to = model_fine$dimension, by = 2 * model_fine$ngrid)
  index1 <- NULL
  index2 <- NULL
  index3 <- NULL
  index4 <- NULL
  for (m in 1:model_coarse$ngrid){
    index1 <- c(index1, seq(from = index_start[m], to = index_start[m]+model_fine$ngrid-2, by = 2))
    index2 <- c(index2, seq(from = index_start[m]+1, to = index_start[m]+model_fine$ngrid-1, by = 2))
    index3 <- c(index3, seq(from = index_start[m]+model_fine$ngrid, to = index_end[m]-1, by = 2))
    index4 <- c(index4, seq(from = index_start[m]+model_fine$ngrid+1, to = index_end[m], by = 2))
  }
  rnorm_coarse <- function(rnorm_fine){
    return((rnorm_fine[index1] + rnorm_fine[index2] + rnorm_fine[index3] + rnorm_fine[index4]) / sqrt(4))
  }

  # multilevel HMC kernel
  multilevel_hmc_kernel <- function(chain_state_coarse, current_pdf_coarse, chain_state_fine, current_pdf_fine){
    current_v_fine <- rnorm(model_fine$dimension) # common velocity or momentum, shared by both multilevel chains
    current_v_coarse <- rnorm_coarse(current_v_fine)
    leapfrog_result_coarse <- leapfrog_integrator(model_coarse, chain_state_coarse, current_v_coarse)
    leapfrog_result_fine <- leapfrog_integrator(model_fine, chain_state_fine, current_v_fine)
    proposed_v_coarse <- - leapfrog_result_coarse$v
    proposed_x_coarse <- leapfrog_result_coarse$x
    proposed_v_fine <- - leapfrog_result_fine$v
    proposed_x_fine <- leapfrog_result_fine$x

    proposed_pdf_coarse <- model_coarse$logtarget(proposed_x_coarse)
    proposed_pdf_fine <- model_fine$logtarget(proposed_x_fine)
    accept_ratio_coarse <- proposed_pdf_coarse - current_pdf_coarse
    accept_ratio_fine <- proposed_pdf_fine - current_pdf_fine
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    kinetic_energy_coarse <- sum(current_v_coarse^2) / 2
    kinetic_energy_fine <- sum(current_v_fine^2) / 2
    accept_ratio_coarse <- accept_ratio_coarse + kinetic_energy_coarse - sum(proposed_v_coarse^2) / 2
    accept_ratio_fine <- accept_ratio_fine + kinetic_energy_fine - sum(proposed_v_fine^2) / 2
    logu <- log(runif(1)) # shared by both chains
    accept_coarse <- FALSE
    accept_fine <- FALSE
    if (is.finite(accept_ratio_coarse)){
      accept_coarse <- (logu < accept_ratio_coarse)
    }

    if (accept_coarse){
      chain_state_coarse <- proposed_x_coarse
      current_pdf_coarse <- proposed_pdf_coarse
    }

    if (is.finite(accept_ratio_fine)){
      accept_fine <- (logu < accept_ratio_fine)
    }

    if (accept_fine){
      chain_state_fine <- proposed_x_fine
      current_pdf_fine <- proposed_pdf_fine
    }

    return(list(chain_state_coarse = chain_state_coarse, current_pdf_coarse = current_pdf_coarse,
                chain_state_fine = chain_state_fine, current_pdf_fine = current_pdf_fine,
                accept_coarse = accept_coarse, accept_fine = accept_fine))
  }

  # coupled4 HMC kernel
  coupled4_hmc_kernel <- function(chain_state_coarse1, chain_state_coarse2,
                                  current_pdf_coarse1, current_pdf_coarse2,
                                  chain_state_fine1, chain_state_fine2,
                                  current_pdf_fine1, current_pdf_fine2){
    current_v_fine <- rnorm(model_fine$dimension) # common velocity or momentum, shared by both multilevel chains
    current_v_coarse <- rnorm_coarse(current_v_fine)

    leapfrog_result_coarse1 <- leapfrog_integrator(model_coarse, chain_state_coarse1, current_v_coarse)
    leapfrog_result_coarse2 <- leapfrog_integrator(model_coarse, chain_state_coarse2, current_v_coarse)
    leapfrog_result_fine1 <- leapfrog_integrator(model_fine, chain_state_fine1, current_v_fine)
    leapfrog_result_fine2 <- leapfrog_integrator(model_fine, chain_state_fine2, current_v_fine)

    proposed_v_coarse1 <- - leapfrog_result_coarse1$v
    proposed_v_coarse2 <- - leapfrog_result_coarse2$v

    proposed_x_coarse1 <- leapfrog_result_coarse1$x
    proposed_x_coarse2 <- leapfrog_result_coarse2$x

    proposed_v_fine1 <- - leapfrog_result_fine1$v
    proposed_v_fine2 <- - leapfrog_result_fine2$v

    proposed_x_fine1 <- leapfrog_result_fine1$x
    proposed_x_fine2 <- leapfrog_result_fine2$x

    proposed_pdf_coarse1 <- model_coarse$logtarget(proposed_x_coarse1)
    proposed_pdf_coarse2 <- model_coarse$logtarget(proposed_x_coarse2)

    proposed_pdf_fine1 <- model_fine$logtarget(proposed_x_fine1)
    proposed_pdf_fine2 <- model_fine$logtarget(proposed_x_fine2)

    accept_ratio_coarse1 <- proposed_pdf_coarse1 - current_pdf_coarse1
    accept_ratio_coarse2 <- proposed_pdf_coarse2 - current_pdf_coarse2

    accept_ratio_fine1 <- proposed_pdf_fine1 - current_pdf_fine1
    accept_ratio_fine2 <- proposed_pdf_fine2 - current_pdf_fine2

    # the acceptance ratio also features the "kinetic energy" term of the extended target
    kinetic_energy_coarse <- sum(current_v_coarse^2) / 2
    kinetic_energy_fine <- sum(current_v_fine^2) / 2

    accept_ratio_coarse1 <- accept_ratio_coarse1 + kinetic_energy_coarse - sum(proposed_v_coarse1^2) / 2
    accept_ratio_coarse2 <- accept_ratio_coarse2 + kinetic_energy_coarse - sum(proposed_v_coarse2^2) / 2

    accept_ratio_fine1 <- accept_ratio_fine1 + kinetic_energy_fine - sum(proposed_v_fine1^2) / 2
    accept_ratio_fine2 <- accept_ratio_fine2 + kinetic_energy_fine - sum(proposed_v_fine2^2) / 2

    logu <- log(runif(1)) # shared by all four chains
    accept_coarse1 <- FALSE
    accept_coarse2 <- FALSE
    accept_fine1 <- FALSE
    accept_fine2 <- FALSE

    if (is.finite(accept_ratio_coarse1)){
      accept_coarse1 <- (logu < accept_ratio_coarse1)
    }
    if (is.finite(accept_ratio_coarse2)){
      accept_coarse2 <- (logu < accept_ratio_coarse2)
    }
    if (is.finite(accept_ratio_fine1)){
      accept_fine1 <- (logu < accept_ratio_fine1)
    }
    if (is.finite(accept_ratio_fine2)){
      accept_fine2 <- (logu < accept_ratio_fine2)
    }

    if (accept_coarse1){
      chain_state_coarse1 <- proposed_x_coarse1
      current_pdf_coarse1 <- proposed_pdf_coarse1
    }
    if (accept_coarse2){
      chain_state_coarse2 <- proposed_x_coarse2
      current_pdf_coarse2 <- proposed_pdf_coarse2
    }
    if (accept_fine1){
      chain_state_fine1 <- proposed_x_fine1
      current_pdf_fine1 <- proposed_pdf_fine1
    }
    if (accept_fine2){
      chain_state_fine2 <- proposed_x_fine2
      current_pdf_fine2 <- proposed_pdf_fine2
    }

    return(list(chain_state_coarse1 = chain_state_coarse1, chain_state_coarse2 = chain_state_coarse2,
                current_pdf_coarse1 = current_pdf_coarse1, current_pdf_coarse2 = current_pdf_coarse2,
                chain_state_fine1 = chain_state_fine1, chain_state_fine2 = chain_state_fine2,
                current_pdf_fine1 = current_pdf_fine1, current_pdf_fine2 = current_pdf_fine2,
                accept_coarse1 = accept_coarse1, accept_coarse2 = accept_coarse2,
                accept_fine1 = accept_fine1, accept_fine2 = accept_fine2))
  }

  # multilevel RWMH kernel
  multilevel_rwmh_kernel <- function(chain_state_coarse, current_pdf_coarse, chain_state_fine, current_pdf_fine){
    rnorm_fine <- rnorm(model_fine$dimension)
    rnorm_coarse <- rnorm_coarse(rnorm_fine)
    proposal_coarse <- chain_state_coarse + proposal_sd * rnorm_coarse
    proposal_fine <- chain_state_fine + proposal_sd * rnorm_fine

    proposal_pdf_coarse <- model_coarse$logtarget(proposal_coarse)
    proposal_pdf_fine <- model_fine$logtarget(proposal_fine)

    logu <- log(runif(1))
    accept_coarse <- FALSE
    accept_fine <- FALSE
    if (is.finite(proposal_pdf_coarse)){
      accept_coarse <- (logu < (proposal_pdf_coarse - current_pdf_coarse))
    }
    if (is.finite(proposal_pdf_fine)){
      accept_fine <- (logu < (proposal_pdf_fine - current_pdf_fine))
    }
    if (accept_coarse){
      chain_state_coarse <- proposal_coarse
      current_pdf_coarse <- proposal_pdf_coarse
    }
    if (accept_fine){
      chain_state_fine <- proposal_fine
      current_pdf_fine <- proposal_pdf_fine
    }

    return(list(chain_state_coarse = chain_state_coarse, current_pdf_coarse = current_pdf_coarse,
                chain_state_fine = chain_state_fine, current_pdf_fine = current_pdf_fine,
                accept_coarse = accept_coarse, accept_fine = accept_fine))

  }

  reflectionmaximal4_coupling <- function(chain_state_coarse1, chain_state_coarse2,
                                          chain_state_fine1, chain_state_fine2){
    ## sample 2-way maximal coupling for fine chain
    # sample first proposal
    randn_fine1 <- rnorm(model_fine$dimension)
    fine1 <- chain_state_fine1 + proposal_sd * randn_fine1

    # difference
    zdiff_fine <- (chain_state_fine1 - chain_state_fine2) / proposal_sd

    # evaluate proposal transition densities at first proposal
    pdf1 <- sum(dnorm(randn_fine1, log = TRUE))
    pdf2 <- sum(dnorm(randn_fine1 + zdiff_fine, log = TRUE))
    logacceptprob <- min(pdf1, pdf2) - pdf1
    loguniform <- log(runif(1))

    if (loguniform < logacceptprob){
      # return common proposal for both chains
      fine2 <- fine1

    } else {
      # perform reflection for second proposal
      evec_fine <- zdiff_fine / sqrt(sum(zdiff_fine^2))
      randn_fine2 <- randn_fine1 - 2 * sum(evec_fine * randn_fine1) * evec_fine
      fine2 <- chain_state_fine2 + proposal_sd * randn_fine2

    }

    ## sample 2-way maximal coupling for coarse chain
    # sample first proposal
    randn_coarse1 <- rnorm_coarse(randn_fine1)
    coarse1 <- chain_state_coarse1 + proposal_sd * randn_coarse1

    # difference
    zdiff_coarse <- (chain_state_coarse1 - chain_state_coarse2) / proposal_sd

    # evaluate proposal transition densities at first proposal
    pdf1 <- sum(dnorm(randn_coarse1, log = TRUE))
    pdf2 <- sum(dnorm(randn_coarse1 + zdiff_coarse, log = TRUE))
    logacceptprob <- min(pdf1, pdf2) - pdf1

    if (loguniform < logacceptprob){
      # return common proposal for both chains
      coarse2 <- coarse1

    } else {
      # perform reflection for second proposal
      evec_coarse <- zdiff_coarse / sqrt(sum(zdiff_coarse^2))
      randn_coarse2 <- randn_coarse1 - 2 * sum(evec_coarse * randn_coarse1) * evec_coarse
      coarse2 <- chain_state_coarse2 + proposal_sd * randn_coarse2

    }

    return(list(coarse1 = coarse1, coarse2 = coarse2, fine1 = fine1, fine2 = fine2))

  }


  # coupled4 RWMH kernel
  coupled4_rwmh_kernel <- function(chain_state_coarse1, chain_state_coarse2,
                                   current_pdf_coarse1, current_pdf_coarse2,
                                   chain_state_fine1, chain_state_fine2,
                                   current_pdf_fine1, current_pdf_fine2){
    # sample from maximal coupling
    proposal_value <- reflectionmaximal4_coupling(chain_state_coarse1, chain_state_coarse2,
                                                  chain_state_fine1, chain_state_fine2)

    proposal_coarse1 <- proposal_value$coarse1
    proposal_coarse2 <- proposal_value$coarse2
    proposal_fine1 <- proposal_value$fine1
    proposal_fine2 <- proposal_value$fine2

    proposal_pdf_coarse1 <- model_coarse$logtarget(proposal_coarse1)
    proposal_pdf_coarse2 <- model_coarse$logtarget(proposal_coarse2)
    proposal_pdf_fine1 <- model_fine$logtarget(proposal_fine1)
    proposal_pdf_fine2 <- model_fine$logtarget(proposal_fine2)

    logu <- log(runif(1))
    accept_coarse1 <- FALSE
    accept_coarse2 <- FALSE
    accept_fine1 <- FALSE
    accept_fine2 <- FALSE

    if (is.finite(proposal_pdf_coarse1)){
      accept_coarse1 <- (logu < (proposal_pdf_coarse1 - current_pdf_coarse1))
    }
    if (is.finite(proposal_pdf_coarse2)){
      accept_coarse2 <- (logu < (proposal_pdf_coarse2 - current_pdf_coarse2))
    }
    if (is.finite(proposal_pdf_fine1)){
      accept_fine1 <- (logu < (proposal_pdf_fine1 - current_pdf_fine1))
    }
    if (is.finite(proposal_pdf_fine2)){
      accept_fine2 <- (logu < (proposal_pdf_fine2 - current_pdf_fine2))
    }

    if (accept_coarse1){
      chain_state_coarse1 <- proposal_coarse1
      current_pdf_coarse1 <- proposal_pdf_coarse1
    }
    if (accept_coarse2){
      chain_state_coarse2 <- proposal_coarse2
      current_pdf_coarse2 <- proposal_pdf_coarse2
    }
    if (accept_fine1){
      chain_state_fine1 <- proposal_fine1
      current_pdf_fine1 <- proposal_pdf_fine1
    }
    if (accept_fine2){
      chain_state_fine2 <- proposal_fine2
      current_pdf_fine2 <- proposal_pdf_fine2
    }

    return(list(chain_state_coarse1 = chain_state_coarse1, chain_state_coarse2 = chain_state_coarse2,
                current_pdf_coarse1 = current_pdf_coarse1, current_pdf_coarse2 = current_pdf_coarse2,
                chain_state_fine1 = chain_state_fine1, chain_state_fine2 = chain_state_fine2,
                current_pdf_fine1 = current_pdf_fine1, current_pdf_fine2 = current_pdf_fine2,
                accept_coarse1 = accept_coarse1, accept_coarse2 = accept_coarse2,
                accept_fine1 = accept_fine1, accept_fine2 = accept_fine2))
  }


  # mixture of multilevel HMC and RWMH kernels
  multilevel_mixture_kernel <- function(chain_state_coarse, current_pdf_coarse, chain_state_fine, current_pdf_fine){
    if (runif(1) < probability_rwmh){
      return(multilevel_rwmh_kernel(chain_state_coarse, current_pdf_coarse, chain_state_fine, current_pdf_fine))
    } else{
      return(multilevel_hmc_kernel(chain_state_coarse, current_pdf_coarse, chain_state_fine, current_pdf_fine))
    }
  }

  # mixture of coupled2 HMC and RWMH kernels
  mixture_coupled4_kernel <- function(chain_state_coarse1, chain_state_coarse2,
                                      current_pdf_coarse1, current_pdf_coarse2,
                                      chain_state_fine1, chain_state_fine2,
                                      current_pdf_fine1, current_pdf_fine2){
    if (runif(1) < probability_rwmh){
      return(coupled4_rwmh_kernel(chain_state_coarse1, chain_state_coarse2,
                                  current_pdf_coarse1, current_pdf_coarse2,
                                  chain_state_fine1, chain_state_fine2,
                                  current_pdf_fine1, current_pdf_fine2))
    } else {
      return(coupled4_hmc_kernel(chain_state_coarse1, chain_state_coarse2,
                                 current_pdf_coarse1, current_pdf_coarse2,
                                 chain_state_fine1, chain_state_fine2,
                                 current_pdf_fine1, current_pdf_fine2))
    }
  }

  return(list(multilevel_hmc_kernel = multilevel_hmc_kernel, coupled4_hmc_kernel = coupled4_hmc_kernel,
              multilevel_kernel = multilevel_mixture_kernel, coupled4_kernel = mixture_coupled4_kernel,
              rnorm_coarse = rnorm_coarse))
}
