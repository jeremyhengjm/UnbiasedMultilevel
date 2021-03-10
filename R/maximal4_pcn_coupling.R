#' @rdname maximal4_pcn_coupling
#' @title maximal coupling of four pCN proposals
#' @description Generate pCN proposals for four chains via maximal coupling of all four chains
#' @param chain_state_coarse1 a vector with coordinates of the first particle at the coarse level
#' @param chain_state_coarse2 a vector with coordinates of the second particle at the coarse level
#' @param chain_state_fine1 a vector with coordinates of the first particle at the fine level
#' @param chain_state_fine2 a vector with coordinates of the second particle at the fine level
#' @param identical_coarse a flag that is True if chains at the coarse level are identical (coincide) and False otherwise
#' @param identical_fine a flag that is True if chains at the fine level are identical (coincide) and False otherwise
#' @param tuning_coarse a list that contains parameters needed for genertaion of pCN proposal at the coarse level: standard devation and rho
#' @param tuning_coarse a list that contains parameters needed for genertaion of pCN proposal at the fine level: standard devation and rho
#'@return a list that contains state of the first chain at the coarse level, state of the second chain at the coarse level, state of the second chain, updated value of the flag "identical", cost of proposal generation
#'@export

maximal4_pcn_coupling <- function(chain_state_coarse1, chain_state_coarse2,
                                  chain_state_fine1, chain_state_fine2,
                                  identical_coarse, identical_fine,
                                  tuning_coarse, tuning_fine){

  cost = 0  # cost of proposal generation
  # tuning parameters that define autoregressive proposal
  proposal_sd_coarse <- tuning_coarse$proposal_sd
  proposal_rho_coarse <- tuning_coarse$proposal_rho
  proposal_sd_factor_coarse <- sqrt(1-proposal_rho_coarse^2) * proposal_sd_coarse

  proposal_sd_fine <- tuning_fine$proposal_sd
  proposal_rho_fine <- tuning_fine$proposal_rho
  proposal_sd_factor_fine <- sqrt(1-proposal_rho_fine^2) * proposal_sd_fine

  # sample first proposal of coarse chain
  coarse1 <- proposal_rho_coarse * chain_state_coarse1 + proposal_sd_factor_coarse * rnorm(dimension)

  # evaluate proposal transition densities at first proposal
  pdf_coarse1 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
  pdf_coarse2 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
  pdf_fine1 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
  pdf_fine2 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
  logacceptprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_coarse1

  if (log(runif(1)) < logacceptprob){
    # return common proposal for all 4 chains
    return(list(coarse1 = coarse1, coarse2 = coarse1,
                fine1 = coarse1, fine2 = coarse1,
                identical_coarse = TRUE, identical_fine = TRUE, cost = cost))

  } else {
    if (identical_coarse){
      # if coarse chains have met: sample common proposal
      reject <- TRUE
      coarse1 <- NA
      while (reject){
        # sample first proposal of coarse chain
        coarse1 <- proposal_rho_coarse * chain_state_coarse1 + proposal_sd_factor_coarse * rnorm(dimension)

        # evaluate proposal transition densities at first proposal
        pdf_coarse1 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_coarse1

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }

      # set second proposal as first so coarse chains stay faithful
      coarse2 <- coarse1
      identical_coarse <- TRUE


    } else {
      # if coarse chains have not met: sample proposals independently
      identical_coarse <- FALSE

      # sample proposal for first coarse chain
      reject <- TRUE
      coarse1 <- NA
      while (reject){
        # sample first proposal of coarse chain
        coarse1 <- proposal_rho_coarse * chain_state_coarse1 + proposal_sd_factor_coarse * rnorm(dimension)

        # evaluate proposal transition densities at first proposal
        pdf_coarse1 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(coarse1, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(coarse1, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_coarse1

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }

      # sample proposal for second coarse chain
      reject <- TRUE
      coarse2 <- NA
      while (reject){
        # sample second proposal of coarse chain
        coarse2 <- proposal_rho_coarse * chain_state_coarse2 + proposal_sd_factor_coarse * rnorm(dimension)

        # evaluate proposal transition densities at second proposal
        pdf_coarse1 <- sum(dnorm(coarse2, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(coarse2, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(coarse2, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(coarse2, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_coarse2

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }

    }

    if (identical_fine){
      # if fine chains have met: sample common proposal
      reject <- TRUE
      fine1 <- NA
      while (reject){
        # sample first proposal of fine chain
        fine1 <- proposal_rho_fine * chain_state_fine1 + proposal_sd_factor_fine * rnorm(dimension)

        # evaluate proposal transition densities at first proposal
        pdf_coarse1 <- sum(dnorm(fine1, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(fine1, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(fine1, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(fine1, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_fine1

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }

      # set second proposal as first so fine chains stay faithful
      fine2 <- fine1
      identical_fine <- TRUE


    } else {
      # if fine chains have not met: sample proposals independently
      identical_fine <- FALSE

      # sample proposal for first fine chain
      reject <- TRUE
      fine1 <- NA
      while (reject){
        # sample first proposal of fine chain
        fine1 <- proposal_rho_fine * chain_state_fine1 + proposal_sd_factor_fine * rnorm(dimension)

        # evaluate proposal transition densities at first proposal
        pdf_coarse1 <- sum(dnorm(fine1, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(fine1, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(fine1, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(fine1, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_fine1

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }

      # sample proposal for second fine chain
      reject <- TRUE
      fine2 <- NA
      while (reject){
        # sample second proposal of fine chain
        fine2 <- proposal_rho_fine * chain_state_fine2 + proposal_sd_factor_fine * rnorm(dimension)

        # evaluate proposal transition densities at second proposal
        pdf_coarse1 <- sum(dnorm(fine2, mean = proposal_rho_coarse * chain_state_coarse1, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_coarse2 <- sum(dnorm(fine2, mean = proposal_rho_coarse * chain_state_coarse2, sd = proposal_sd_factor_coarse, log = TRUE))
        pdf_fine1 <- sum(dnorm(fine2, mean = proposal_rho_fine * chain_state_fine1, sd = proposal_sd_factor_fine, log = TRUE))
        pdf_fine2 <- sum(dnorm(fine2, mean = proposal_rho_fine * chain_state_fine2, sd = proposal_sd_factor_fine, log = TRUE))
        logrejectprob <- min(pdf_coarse1, pdf_coarse2, pdf_fine1, pdf_fine2) - pdf_fine2

        # accept or reject
        reject <- (log(runif(1)) < logrejectprob)

      }
    }
    return(list(coarse1 = coarse1, coarse2 = coarse2, fine1 = fine1, fine2 = fine2,
                identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))
  }
}

