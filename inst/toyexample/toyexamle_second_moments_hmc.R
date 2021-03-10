# This script can be run in order to generate the data to estimate the rate of decay
# os the second moments of the increments estimated with HMC and RWM reflection maximal coupling
rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)
set.seed(1)
# data generation:
xdata_ref = npyLoad('./toyexample_data.npy')  # you can use only that line of code if you simply want to load the data

su_param = 4.0  # std of prior distribution
sm_param = 1.0  # std of likelihood function

log_lh_fun <- function(u, l, sigma_u = su_param, sigma_m = sm_param,
                       xdata = xdata_ref)
{
    # output: log(likelihood) + log(prior)
    # u - vector of model parameters
    # l - level
    # sigma_u = su_param - std of the prior distribution
    # sigma_m = sm_param - std of the likelihood
    # xdata = xdata_ref - vector of observations
    nu = length(u)
    xmeas = observation_toyexample(u, l)  # vector of observations ~ model predictions ~ measurements
    llh = 0.0  # init variable for output
    # contribution of the prior distribution
    for (k in 1 : nu)
    {
        llh = llh - 0.5 * log(2.0 * 3.1415926) - log(sigma_u)
        llh = llh - 0.5 * u[k] * u[k] / sigma_u / sigma_u
    }

    # contribution of the likelihood
    for (k in 1 : length(xmeas))
    {
        llh = llh - 0.5 * log(2.0 * 3.1415926) - log(sigma_m)
        llh = llh - 0.5 * (xmeas[k] - xdata[k]) * (xmeas[k] - xdata[k]) / sigma_m / sigma_m
    }
    return(llh)
}


grad_log_lh_fun <- function(u, l, sigma_u = su_param, sigma_m = sm_param,
                            xdata = xdata_ref)
{
    # output: derivatives of (log(likelihood) + log(prior))
    # with respect to model parameters: u

    # u - vector of model parameters
    # l - level
    # sigma_u = su_param - std of the prior distribution
    # sigma_m = sm_param - std of the likelihood
    # xdata = xdata_ref - vector of observations

    nu = length(u)
    xmeas = observation_toyexample(u, l)  # vector of observations ~ model predictions ~ measurements
    xmeas_grad = observation_grad_toyexample(u, l)  # derivatives of vector of observations with respect to u
    grad_llh = 0.0 * u  # init variable for output
    # contribution of the prior distribution
    for (k in 1 : nu)
    {
        grad_llh[k] = grad_llh[k] - u[k] / sigma_u / sigma_u
    }
    # contribution of the likelihood
    for (k0 in 1 : nu)
    {
        for (k in 1 : length(xmeas))
        {
            grad_llh[k0] = grad_llh[k0] - (xmeas[k] - xdata[k]) / sigma_m / sigma_m * xmeas_grad[k, k0]
        }
    }
    return(grad_llh)
}


objective_function1 <- function(l, u)
{
    # test function
    return(u)
}

objective_function <- function(l, u, sigma_u = su_param, sigma_m = sm_param,
                               xdata = xdata_ref)
{
    # gradient of log(Bayesian_evidence_factor) with respect to std the likelihood function
    # output: dlog(Bayesian_evidence_factor) / d sigma_m

    # u - vector of model parameters
    # l - level
    # sigma_u = su_param - std of the prior distribution
    # sigma_m = sm_param - std of the likelihood
    # xdata = xdata_ref - vector of observations

    nu = length(u)
    xmeas = observation_toyexample(u, l)  # value of observations for model paramets u and level l
    nv = 3  # dimension of the output vector
    llh_grad <- array(rep(0, nv), dim = c(nv))  # init output

    # contribution of the likelihood
    llh_grad_term = 0.0
    for (k in 1 : length(xmeas))
    {
        llh_grad_term = llh_grad_term + (xmeas[k] - xdata[k]) * (xmeas[k] - xdata[k]) / sigma_m / sigma_m / length(xmeas)
    }
    llh_grad_term = llh_grad_term - 1.0
    llh_grad_term = llh_grad_term * length(xmeas) / sigma_m
    llh_grad[1] = llh_grad[1] + llh_grad_term
    llh_grad[2] = u[1]
    llh_grad[3] = u[2]
    return(llh_grad)
}

init_pl <- function(beta, lmax=20)
{
    # beta = (4 + 1) / 2
    # probability distribution over levels: w[l] = P_L(l - 1)
    # lmax - trunction of the P_L
    # beta - rate of the decay
    w = 1 : lmax
    w = 2.0 ^ (-beta * w)
    w = w / sum(w)
    return(w)
}

get_l_sample <- function(w)
{
    # sample l from P_l
    u = runif(1, min=0.0, max=1.0)
    lmax = length(w)
    l0 = 0
    w0 = 1.0
    wsum = 0.0
    for(l in (1:lmax))
    {
        wsum = sum(w[1 : l])
        if(wsum < u)
        {
            l0 = l
        }
    }
    return(l0)
}

# specify sequence of target distributions
dimension <- 2
logtarget <- function(l, x) log_lh_fun(x, l)
gradlogtarget <- function(l, x) grad_log_lh_fun(x, l)
level <- 0 # level

# vary levels
nrepeats <- 100
constant <- 0.01  # fix tuning parameters across levels for simplicity (could vary this!)
stepsize <- 0.001
nsteps <- 1 + 2 * floor(1 / stepsize) # set integration time to approximately one
probability_maximal_coupling <- 0.10 # probability of selecting the maximal coupling
tuning_coarse <- list(proposal_sd = 1.0e-4, stepsize = 1.0e-1, nsteps = 10)
tuning_fine <- list(proposal_sd = 1.0e-4, stepsize = 1.0e-1, nsteps = 10)
tuning <- list(proposal_sd = 1.0e-4, stepsize = 1.0e-1, nsteps = 10)

rinit <- function(level){
  chain_state <- rnorm(dimension, mean=0.0, sd=su_param)
  current_pdf <- logtarget(level, chain_state)
  chain_data <- list(chain_state = chain_state, current_pdf = current_pdf)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

beta = (4.0 + 1.0) / 2.0  # value of P_L(l) decay
w = init_pl(beta)  # values of P_L(l)

k0 = 100  # k from the paper
m0 = 10 * k0  # m form the paper
niter = 1000  # number os SGD steps
nreal = 1  # number of realisations
nv = 3  # dimension of the paramet space (of the gradient in sgd)
nlevel = 10
gf = 0.2  # factor for the SGD step-size. SGD step size is gf / sgd_step
xval = array(rep(0, nv, dim = c(nv)))

for(kreal in 1 : nreal)
{
    # memory allocation
    sgd_vec <- array(rep(0, niter * nv), dim = c(niter, nv))
    xval = array(rep(0, nv), dim = c(nv))
    nopr <- array(rep(0, niter), dim = c(niter)) # total operations counter
    su_param = 4.0  # std of prior distribution
    sm_param = 1.0  # std of likelihood function

    for(vlevel in 0 : nlevel)
    {
        sgd_vec = 0.0 * sgd_vec
        for(kiter in 1 : niter)
        {
            xval = 0.0 * xval
            if(vlevel == 0)
            {
                x_expect = unbiased_expectation(level = vlevel,
                                                rinit = rinit,
                                                single_kernel = single_hmc_kernel,
                                                tuning = tuning,
                                                coupled_kernel = mixture2_hmc_kernel,
                                                proposal_coupling = reflectionmaximal2_rwmh_coupling,
                                                h = function(l, x) objective_function(l, x),
                                                k = k0, m = m0,
                                                max_iterations = Inf)
                xval <- xval + x_expect$uestimator # / w[1 + vlevel]
                nopr[kiter] = x_expect$cost
            }
            if(vlevel > 0)
            {
                x_increm = unbiased_increment(level = 0 + vlevel,
                                              rinit = rinit,
                                              single_kernel = single_pcn_kernel,
                                              coupled2_kernel = mixture2_hmc_kernel_extra_level,
                                              coupled4_kernel = mixture4_hmc_kernel,
                                              proposal_coupling2 = reflectionmaximal2_rwmh_coupling,
                                              proposal_coupling4 = reflectionmaximal4_rwmh_synchronous_coupling,
                                              tuning = tuning,
                                              tuning_coarse = tuning,
                                              tuning_fine = tuning,
                                              h = function(l, x) objective_function(l, x),
                                              k = k0, m = m0,
                                              sampling_factor = min(0.5, 1.0 / (2.0 ^ (2 * vlevel + 1))),
                                              max_iterations = Inf)
                xval <- xval + x_increm$uestimator # / w[1 + vlevel]
                nopr[kiter] = x_increm$cost
            }
            sgd_vec[kiter, 1 : nv] = xval[1 : nv]
            fname = paste0("./toyexample_second_moment_hmc_level_", vlevel, ".rds")
            saveRDS(sgd_vec, file = fname)
            fname = paste0("./toyexample_second_moment_hmc_level_", vlevel, ".npy")
            npySave(fname, sgd_vec)
        }
    }
}














