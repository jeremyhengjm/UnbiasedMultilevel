#' @rdname synchronous4_pcn_coupling
#' @title Syncronous coupling of four pCN proposals
#' @description Sampling of common increment for all four particles
#' @param chain_state_coarse1 a vector with the position of the first particle at the coarse level (level - 1)
#' @param chain_state_coarse2 a vector with the position of the second particle at the coarse level (level - 1)
#' @param chain_state_fine1 a vector with the position of the first particle at the fine level (level)
#' @param chain_state_fine2 a vector with the position of the second particle at the fine level (level)
#' @param identical_coarse a boolean variable takes the value True if two chains at the coarse level coincide and False otherwise
#' @param identical_fine a boolean variable takes the value True if two chains at the fine level coincide and False otherwise
#' @param tuning_coarse a list that contains parameters required for pCN proposal generation at the coarse level: standard deviation and rho
#' @param tuning_fine a list that contains parameters required for pCN proposal generation at the fine level: standard deviation and rho
#'@return a list with updated values of states of the first and the second chain at the coarse level, updated values of states of the first and the second chain at the fine level, updated value of the "identical_coarse" flag, updated value of the "identical_fine" flag and cost of computations
#'@export

synchronous4_pcn_coupling <- function(chain_state_coarse1, chain_state_coarse2,
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

  # sample common Gaussian increment
  noise <- rnorm(dimension)

  # sample proposals of coarse chains
  coarse1 <- proposal_rho_coarse * chain_state_coarse1 + proposal_sd_factor_coarse * noise
  coarse2 <- proposal_rho_coarse * chain_state_coarse2 + proposal_sd_factor_coarse * noise

  # sample proposals of fine chains
  fine1 <- proposal_rho_fine * chain_state_fine1 + proposal_sd_factor_fine * noise
  fine2 <- proposal_rho_fine * chain_state_fine2 + proposal_sd_factor_fine * noise

  # if chains are identical, they stay together;
  # if chains are not identical, the synchronous coupling cannot allow them to be

  return(list(coarse1 = coarse1, coarse2 = coarse2, fine1 = fine1, fine2 = fine2,
              identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))

}

