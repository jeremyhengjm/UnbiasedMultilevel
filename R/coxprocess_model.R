#' @rdname coxprocess_model
#' @title Construct log-Gaussian Cox process model
#' @description Construct log-Gaussian Cox process model at a given spatial resolution
#' @param parameters a list of model parameters (sigmasq, mu, beta)
#' @param dataset a list of observed spatial locations
#' @param level an integer that determines the spatial resolution
#'@return a list containing the keys \code{dimension}, \code{ngrid}, \code{logtarget}, \code{gradlogtarget}, \code{rinit} and \code{gradient_parameters}
#'@export

coxprocess_model <- function(parameters, dataset, level){

  # spatial resolution
  ngrid <- 2^level
  dimension <- ngrid^2
  grid <- seq(from = 0, to = 1, length.out = ngrid+1)
  area <- 1 / dimension # area of each grid cell

  # compute observed counts in each grid cell
  data_counts <- rep(0, dimension)
  for (i in 1:ngrid){
    for (j in 1:ngrid){
      logical_x <- (dataset$x > grid[i]) * (dataset$x < grid[i+1])
      logical_y <- (dataset$y > grid[j]) * (dataset$y < grid[j+1])
      data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
    }
  }

  # prior distribution
  prior_mean <- rep(parameters$mu, dimension)
  prior_cov <- matrix(nrow = dimension, ncol = dimension)
  for (m in 1:dimension){
    for (n in 1:dimension){
      index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
      index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
      prior_cov[m,n] <- parameters$sigmasq * exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameters$beta) )
    }
  }
  prior_precision <- solve(prior_cov)
  prior_precision_chol <- t(chol(prior_precision))
  prior_cov_chol <- solve(t(prior_precision_chol))
  prior <- list()
  prior$logdensity <- function(x){
    return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
  }
  prior$gradlogdensity <- function(x){
    return(as.numeric((prior_mean - x) %*% prior_precision))
  }

  # likelihood function
  likelihood <- list()
  likelihood$log <- function(x){
    return(coxprocess_loglikelihood(matrix(x, nrow = 1), data_counts, area))
  }
  likelihood$gradlog <- function(x){
    return( data_counts - area * exp(x) )
  }

  # posterior distribution
  logtarget <- function(x) prior$logdensity(x) + likelihood$log(x)
  gradlogtarget <- function(x) prior$gradlogdensity(x) + likelihood$gradlog(x)

  # initial distribution
  # rinit <- function() as.numeric(fast_rmvnorm(1, prior_mean, prior_cov))
  rinit <- function(randn) prior_mean + as.numeric(randn %*% t(prior_cov_chol))

  # functional for stochastic gradient algorithm
  nparameters <- 3
  gradient_parameters <- function(x){
    # precompute
    vector_matrix_product <- (x - parameters$mu) %*% prior_precision

    # preallocate
    gradient <- rep(0, nparameters)
    for (m in 1:dimension){
      for (n in 1:dimension){
        # precompute
        index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
        index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
        norm_index <- sqrt(sum((index_m - index_n)^2))

        # with respect to mu
        gradient[1] <- gradient[1] + (0.5 * x[n] + 0.5 * x[m] - parameters$mu) * prior_precision[m, n]

        # with respect to sigmasq
        dsigmasq <- exp(- norm_index / (ngrid * parameters$beta) )
        gradient[2] <- gradient[2] + 0.5 * vector_matrix_product[n] * vector_matrix_product[m] * dsigmasq -
          0.5 * prior_precision[m,n] * dsigmasq

        # with respect to beta
        dbeta <- parameters$sigmasq * exp(- norm_index / (ngrid * parameters$beta) ) * norm_index / (ngrid * parameters$beta^2)
        gradient[3] <- gradient[3] + 0.5 * vector_matrix_product[n] * vector_matrix_product[m] * dbeta -
          0.5 * prior_precision[m,n] * dbeta
      }
    }
    return(gradient)
  }

  return(list(dimension = dimension, ngrid = ngrid,
              logtarget = logtarget, gradlogtarget = gradlogtarget,
              rinit = rinit, gradient_parameters = gradient_parameters))

}
