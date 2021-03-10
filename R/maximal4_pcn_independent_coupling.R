#' @rdname maximal4_pcn_independent_coupling
#' @title independent maximal coupling of 4 pCN proposals
#' @description Generate pCN proposals for two chains via maximal coupling for pairs of chains at the coarse and fine levels independently
#' @param chain_state_coarse1 a vector with coordinates of the first particle at the coarse level
#' @param chain_state_coarse2 a vector with coordinates of the second particle at the coarse level
#' @param chain_state_fine1 a vector with coordinates of the first particle at the fine level
#' @param chain_state_fine2 a vector with coordinates of the second particle at the fine level
#' @param identical_coarse a flag that is True if chains at the coarse level are identical (coincide) and False otherwise
#' @param identical_fine a flag that is True if chains at the fine level are identical (coincide) and False otherwise
#' @param tuning_coarse a list that contains parameters needed for genertaion of pCN proposal at the coarse level: standard devation and rho
#' @param tuning_coarse a list that contains parameters needed for genertaion of pCN proposal at the fine level: standard devation and rho
#'@return a list that contains state of the first chain at the coarse level, state of the second chain at the coarse level, state of the second chain, updated values of flags "identical_coarse" and "identical_fine", cost of proposal generation
#'@export

maximal4_pcn_independent_coupling <- function(chain_state_coarse1, chain_state_coarse2,
                                              chain_state_fine1, chain_state_fine2,
                                              identical_coarse, identical_fine,
                                              tuning_coarse, tuning_fine){
  cost = 0  # cost of proposal generation
  # sample 2-way maximal coupling for coarse chain
  maximal_output <- maximal2_pcn_coupling(chain_state_coarse1, chain_state_coarse2, identical_coarse, tuning_coarse)
  coarse1 <- maximal_output$state1
  coarse2 <- maximal_output$state2
  identical_coarse <- maximal_output$identical
  cost = cost + maximal_output$cost
  
  # sample 2-way maximal coupling for fine chain
  maximal_output <- maximal2_pcn_coupling(chain_state_fine1, chain_state_fine2, identical_fine, tuning_fine)
  fine1 <- maximal_output$state1
  fine2 <- maximal_output$state2
  identical_fine <- maximal_output$identical
  cost = cost + maximal_output$cost
  
  return(list(coarse1 = coarse1, coarse2 = coarse2, fine1 = fine1, fine2 = fine2,
              identical_coarse = identical_coarse, identical_fine = identical_fine, cost = cost))
  
}

