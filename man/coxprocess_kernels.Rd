% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxprocess_kernels.R
\name{coxprocess_kernels}
\alias{coxprocess_kernels}
\title{Markov chain Monte Carlo kernels for log-Gaussian Cox process model}
\usage{
coxprocess_kernels(model, tuning)
}
\arguments{
\item{model}{list containing \code{logtarget}, \code{gradlogtarget} and \code{dimension}}

\item{tuning}{list containing \code{stepsize} and \code{nsteps}}
}
\value{
a list containing the keys \code{single_hmc_kernel}, \code{coupled2_hmc_kernel}, \code{single_kernel} and \code{coupled2_kernel}
}
\description{
This function defines Markov chain Monte Carlo kernels and its couplings
}
