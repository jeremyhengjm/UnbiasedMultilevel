% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_pcn_kernel.R
\name{single_pcn_kernel}
\alias{single_pcn_kernel}
\title{pCN kernel for a single MCMC chain that target probability distribution at the given level}
\usage{
single_pcn_kernel(level, state, tuning)
}
\arguments{
\item{level}{a integer that determines density of probability distribution inthe multi-level approach}

\item{state}{a list with current position of the particle (element of the chain) and log of the density of probability distribution}

\item{tuning}{a list of parameters for pCN iteration: standard deviation and rho}
}
\value{
a updated state of the chain and cost of computations
}
\description{
Generation of pCN proposal and accept/reject step
}
