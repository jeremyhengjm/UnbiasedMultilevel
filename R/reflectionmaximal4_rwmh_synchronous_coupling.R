#' @rdname reflectionmaximal4_rwmh_synchronous_coupling
#' @title Reflection maximal coupling of four RWM proposals
#' @description Generate RWM proposals for four chains via reflection maximal coupling of all four chains
#' @param chain_state_coarse1 a vector with coordinates of the first particle at the coarse level
#' @param chain_state_coarse2 a vector with coordinates of the second particle at the coarse level
#' @param chain_state_fine1 a vector with coordinates of the first particle at the fine level
#' @param chain_state_fine2 a vector with coordinates of the second particle at the fine level
#' @param identical_coarse a flag that is True if chains at the coarse level are identical (coincide) and False otherwise
#' @param identical_fine a flag that is True if chains at the fine level are identical (coincide) and False otherwise
#' @param tuning_coarse a list that contains parameters needed for genertaion of RWM proposal at the coarse level: standard devation
#' @param tuning_coarse a list that contains parameters needed for genertaion of RWm proposal at the fine level: standard devation
#'@return a list that contains state of the first chain at the coarse level, state of the second chain at the coarse level, state of the second chain, updated value of the flag "identical", cost of proposal generation
#'@export

reflectionmaximal4_rwmh_synchronous_coupling <- function(chain_state_coarse1, chain_state_coarse2,
                                                         chain_state_fine1, chain_state_fine2,
                                                         identical_coarse, identical_fine,
                                                         tuning_coarse, tuning_fine){
  cost = 0  # cost of proposal generation
  # extract tuning parameters
  proposal_sd_coarse <- tuning_coarse$proposal_sd
  proposal_sd_fine <- tuning_fine$proposal_sd
  
  ## sample 2-way maximal coupling for coarse chain
  # sample first proposal
  randn_coarse1 <- rnorm(dimension)
  coarse1 <- chain_state_coarse1 + proposal_sd_coarse * randn_coarse1
  
  # difference
  zdiff_coarse <- (chain_state_coarse1 - chain_state_coarse2) / proposal_sd_coarse
  
  # evaluate proposal transition densities at first proposal
  pdf1 <- sum(dnorm(randn_coarse1, log = TRUE))
  pdf2 <- sum(dnorm(randn_coarse1 + zdiff_coarse, log = TRUE))
  logacceptprob <- min(pdf1, pdf2) - pdf1
  loguniform <- log(runif(1))
  
  if (loguniform < logacceptprob){
    # return common proposal for both chains
    coarse2 <- coarse1
    identical_coarse <- TRUE
    
  } else {
    # perform reflection for second proposal
    evec_coarse <- zdiff_coarse / sqrt(sum(zdiff_coarse^2))
    randn_coarse2 <- randn_coarse1 - 2 * sum(evec_coarse * randn_coarse1) * evec_coarse
    coarse2 <- chain_state_coarse2 + proposal_sd_coarse * randn_coarse2
    identical_coarse <- FALSE
    
  }
  
  ## sample 2-way maximal coupling for fine chain
  # sample first proposal
  randn_fine1 <- randn_coarse1 # synchronous coupling
  fine1 <- chain_state_fine1 + proposal_sd_fine * randn_fine1
  
  # difference
  zdiff_fine <- (chain_state_fine1 - chain_state_fine2) / proposal_sd_fine
  
  # evaluate proposal transition densities at first proposal
  pdf2 <- sum(dnorm(randn_fine1 + zdiff_fine, log = TRUE))
  logacceptprob <- min(pdf1, pdf2) - pdf1
  
  if (loguniform < logacceptprob){
    # return common proposal for both chains
    fine2 <- fine1
    identical_fine <- TRUE
    
  } else {
    # perform reflection for second proposal
    evec_fine <- zdiff_fine / sqrt(sum(zdiff_fine^2))
    randn_fine2 <- randn_fine1 - 2 * sum(evec_fine * randn_fine1) * evec_fine
    fine2 <- chain_state_fine2 + proposal_sd_fine * randn_fine2
    identical_fine <- FALSE
    
  }
  
  return(list(coarse1 = coarse1, coarse2 = coarse2, fine1 = fine1, fine2 = fine2,
              identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))
  
}

