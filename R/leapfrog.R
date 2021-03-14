#' @rdname leapfrog
#' @title Leapfrog inregration for HMC kernel
#' @description leapfrog integrtion
#' @param level an integer that controls probability distribution utilized for calculation of gradients of the potential
#' @param x position of the particle
#' @param v velocity of the particle
#' @param tuning a list that contains the stepsize and number of steps for the leapfrog integration
#'@return list that contains new position of the particle, new velocsity of the particle and the cost of leapfrog integration (proportional to number of gradient evaluation)
#'@export

leapfrog <- function(level, x, v, tuning){

  # extract tuning parameters
  stepsize <- tuning$stepsize
  nsteps <- tuning$nsteps
  cost <- 0
  dimension <- length(v)

  v <- v + stepsize * gradlogtarget(level, x) / 2
  cost <- cost + dimension * 2^level
  for (step in 1:nsteps){
    x <- x + stepsize * v
    if (step != nsteps){
      v <- v + stepsize * gradlogtarget(level, x)
      cost <- cost + dimension * 2^level
    }
  }
  v <- v + stepsize * gradlogtarget(level, x) / 2
  cost <- cost + length(v)*2^level

  # we could negate the momentum but we don't use it here
  return(list(x = x, v = v, cost = cost))
}
