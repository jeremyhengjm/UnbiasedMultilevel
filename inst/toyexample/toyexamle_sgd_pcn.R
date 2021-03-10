# This script can be run in order to generate the data to perform MAP estimate of the hyperparameter
# reflection maximal coupling of pCN proposals is utilized to estimate the gradient

rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)
set.seed(2)

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
    ny = length(u)
    xmeas = observation_toyexample(u, l)  # value of observations for model paramets u and level l
    nv = 1  # dimension of the output vector
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
proposal_sd = 0.0 + su_param
proposal_rho = 0.95
tuning <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)
tuning_fine <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)
tuning_coarse <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)

rinit <- function(level){
  chain_state <- rnorm(dimension, mean=0.0, sd=su_param)
  current_pdf <- logtarget(level, chain_state)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

beta = (4.0 + 1.0) / 2.0  # value of P_L(l) decay
w = init_pl(beta)  # values of P_L(l)

k0 = 200  # k from the paper
m0 = 2 * k0  # m form the paper
niter = 1000  # number os SGD steps
nreal = 100  # number of realisations
nv = 1  # dimension of the paramet space (of the gradient in sgd)
nlevel = 10

gf_case = 1

if(gf_case == 1)
{
    gf = 0.03  # factor for the SGD step-size. SGD step size is gf / sgd_step
    xval = array(rep(0, nv, dim = c(nv)))

    for(kreal in 1 : nreal)
    {
        # memory allocation
        sgd_vec <- array(rep(0, niter * nv), dim = c(niter, nv))
        sigma_rep <- array(rep(0, niter), dim = c(niter))  # array of sigma values
        xval = array(rep(0, nv), dim = c(nv))
        nopr <- array(rep(0, niter), dim = c(niter)) # total operations counter
        su_param = 4.0  # std of prior distribution
        sm_param = 1.0  # std of likelihood function

        for(kiter in 1 : niter)
        {
            sgd_vec = 0.0 * sgd_vec
            xval = 0.0 * xval
            vlevel = get_l_sample(w)
            if(vlevel == 0)
            {
                x_expect = unbiased_expectation(level = vlevel,
                                                rinit = rinit,
                                                single_kernel = single_pcn_kernel,
                                                tuning = tuning,
                                                coupled_kernel = coupled2_pcn_kernel,
                                                proposal_coupling = reflectionmaximal2_pcn_coupling,
                                                h = function(l, x) objective_function(l, x),
                                                k = k0, m = m0,
                                                max_iterations = Inf)
                xval <- xval + x_expect$uestimator / w[1 + vlevel]
                nopr[kiter] = x_expect$cost
            }
            if(vlevel > 0)
            {
                x_increm = unbiased_increment(level = 0 + vlevel,
                                              rinit = rinit,
                                              single_kernel = single_pcn_kernel,
                                              coupled2_kernel = coupled2_pcn_kernel_extra_level,
                                              coupled4_kernel = coupled4_pcn_kernel,
                                              proposal_coupling2 = synchronous2_pcn_coupling,
                                              proposal_coupling4 = reflectionmaximal4_pcn_synchronous_coupling,
                                              tuning = tuning,
                                              tuning_coarse = tuning,
                                              tuning_fine = tuning,
                                              h = function(l, x) objective_function(l, x),
                                              k = k0, m = m0,
                                              sampling_factor = min(0.5, 1.0 / (2.0 ^ (2 * vlevel + 1))),
                                              max_iterations = Inf)
                xval <- xval + x_increm$uestimator / w[1 + vlevel]
                nopr[kiter] = x_increm$cost
            }
            sgd_vec[kiter, 1 : nv] = xval[1 : nv]
            # regularization for to secure from the blow-out at the early steps
            if((abs(gf * sgd_vec[kiter, 1]) * sm_param) > 3.0)
            {
                #print("high gradient correction")
                sgd_vec[kiter, 1] = 3.0 / sm_param / gf * sgd_vec[kiter, 1] / abs(sgd_vec[kiter, 1])
            }
            # update the std in the likelihood function
            sm_param = sm_param * exp(gf / kiter * sgd_vec[kiter, 1] * sm_param)

            # record the data
            sigma_rep[kiter] = 0.0 + sm_param

            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v7.rds")
            saveRDS(sigma_rep, file = fname)
            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v7.npy")
            npySave(fname, sigma_rep)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v7.rds")
            saveRDS(sgd_vec, file = fname)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v7.npy")
            npySave(fname, sgd_vec)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v7.rds")
            saveRDS(nopr, file = fname)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v7.npy")
            npySave(fname, nopr)

        }
    }
}



if(gf_case == 2)
{
    gf = 0.01  # factor for the SGD step-size. SGD step size is gf / sgd_step
    xval = array(rep(0, nv, dim = c(nv)))

    for(kreal in 1 : nreal)
    {
        # memory allocation
        sgd_vec <- array(rep(0, niter * nv), dim = c(niter, nv))
        sigma_rep <- array(rep(0, niter), dim = c(niter))  # array of sigma values
        xval = array(rep(0, nv), dim = c(nv))
        nopr <- array(rep(0, niter), dim = c(niter)) # total operations counter
        su_param = 4.0  # std of prior distribution
        sm_param = 1.0  # std of likelihood function

        for(kiter in 1 : niter)
        {
            sgd_vec = 0.0 * sgd_vec
            xval = 0.0 * xval
            writeLines('0', "nopr.txt")  # set counter to zero
            vlevel = get_l_sample(w)
            if(vlevel == 0)
            {
                x_expect = unbiased_expectation(level = vlevel,
                                                rinit = rinit,
                                                single_kernel = single_pcn_kernel,
                                                tuning = tuning,
                                                coupled_kernel = coupled2_pcn_kernel,
                                                proposal_coupling = reflectionmaximal2_pcn_coupling,
                                                h = function(l, x) objective_function(l, x),
                                                k = k0, m = m0,
                                                max_iterations = Inf)
                xval <- xval + x_expect$uestimator / w[1 + vlevel]
                nopr[kiter] = x_expect$cost
            }
            if(vlevel > 0)
            {
                x_increm = unbiased_increment(level = 0 + vlevel,
                                              rinit = rinit,
                                              single_kernel = single_pcn_kernel,
                                              coupled2_kernel = coupled2_pcn_kernel_extra_level,
                                              coupled4_kernel = coupled4_pcn_kernel,
                                              proposal_coupling2 = synchronous2_pcn_coupling,
                                              proposal_coupling4 = reflectionmaximal4_pcn_synchronous_coupling,
                                              tuning = tuning,
                                              tuning_coarse = tuning,
                                              tuning_fine = tuning,
                                              h = function(l, x) objective_function(l, x),
                                              k = k0, m = m0,
                                              sampling_factor = min(0.5, 1.0 / (2.0 ^ (2 * vlevel + 1))),
                                              max_iterations = Inf)
                xval <- xval + x_increm$uestimator / w[1 + vlevel]
                nopr[kiter] = x_expect$cost
            }
            sgd_vec[kiter, 1 : nv] = xval[1 : nv]
            # regularization for to secure from the blow-out at the early steps
            if((abs(gf * sgd_vec[kiter, 1]) * sm_param) > 1.0)
            {
                #print("high gradient correction")
                sgd_vec[kiter, 1] = 1.0 / sm_param / gf * sgd_vec[kiter, 1] / abs(sgd_vec[kiter, 1])
            }
            # update the std in the likelihood function
            sm_param = sm_param * exp(gf / kiter * sgd_vec[kiter, 1] * sm_param)

            # record the data
            sigma_rep[kiter] = 0.0 + sm_param

            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v2.rds")
            saveRDS(sigma_rep, file = fname)
            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v2.npy")
            npySave(fname, sigma_rep)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v2.rds")
            saveRDS(sgd_vec, file = fname)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v2.npy")
            npySave(fname, sgd_vec)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v2.rds")
            saveRDS(nopr, file = fname)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v2.npy")
            npySave(fname, nopr)
        }
    }
}


if(gf_case == 3)
{
    gf = 0.003 # factor for the SGD step-size. SGD step size is gf / sgd_step
    xval = array(rep(0, nv, dim = c(nv)))

    for(kreal in 1 : nreal)
    {
        # memory allocation
        sgd_vec <- array(rep(0, niter * nv), dim = c(niter, nv))
        sigma_rep <- array(rep(0, niter), dim = c(niter))  # array of sigma values
        xval = array(rep(0, nv), dim = c(nv))
        nopr <- array(rep(0, niter), dim = c(niter)) # total operations counter
        su_param = 4.0  # std of prior distribution
        sm_param = 1.0  # std of likelihood function

        for(kiter in 1 : niter)
        {
            sgd_vec = 0.0 * sgd_vec
            xval = 0.0 * xval
            vlevel = get_l_sample(w)
            if(vlevel == 0)
            {
                x_expect = unbiased_expectation(level = vlevel,
                                                rinit = rinit,
                                                single_kernel = single_pcn_kernel,
                                                tuning = tuning,
                                                coupled_kernel = coupled2_pcn_kernel,
                                                proposal_coupling = reflectionmaximal2_pcn_coupling,
                                                h = function(l, x) objective_function(l, x),
                                                k = k0, m = m0,
                                                max_iterations = Inf)
                xval <- xval + x_expect$uestimator / w[1 + vlevel]
                nopr[kiter] = x_expect$cost
            }
            if(vlevel > 0)
            {
                x_increm = unbiased_increment(level = 0 + vlevel,
                                              rinit = rinit,
                                              single_kernel = single_pcn_kernel,
                                              coupled2_kernel = coupled2_pcn_kernel_extra_level,
                                              coupled4_kernel = coupled4_pcn_kernel,
                                              proposal_coupling2 = synchronous2_pcn_coupling,
                                              proposal_coupling4 = reflectionmaximal4_pcn_synchronous_coupling,
                                              tuning = tuning,
                                              tuning_coarse = tuning,
                                              tuning_fine = tuning,
                                              h = function(l, x) objective_function(l, x),
                                              k = k0, m = m0,
                                              sampling_factor = min(0.5, 1.0 / (2.0 ^ (2 * vlevel + 1))),
                                              max_iterations = Inf)
                xval <- xval + x_increm$uestimator / w[1 + vlevel]
                nopr[kiter] = x_incremt$cost
            }
            sgd_vec[kiter, 1 : nv] = xval[1 : nv]
            # regularization for to secure from the blow-out at the early steps
            if((abs(gf * sgd_vec[kiter, 1]) * sm_param) > 1.0)
            {
                #print("high gradient correction")
                sgd_vec[kiter, 1] = 1.0 / sm_param / gf * sgd_vec[kiter, 1] / abs(sgd_vec[kiter, 1])
            }
            # update the std in the likelihood function
            sm_param = sm_param * exp(gf / kiter * sgd_vec[kiter, 1] * sm_param)

            # record the data
            sigma_rep[kiter] = 0.0 + sm_param

            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v6.rds")
            saveRDS(sigma_rep, file = fname)
            fname = paste0("./toyexample_sgd_sigma_pcn_kreal_", kreal, "_v6.npy")
            npySave(fname, sigma_rep)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v6.rds")
            saveRDS(sgd_vec, file = fname)
            fname = paste0("./toyexample_sgd_pcn_nreal_", kreal, "_v6.npy")
            npySave(fname, sgd_vec)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v6.rds")
            saveRDS(nopr, file = fname)
            fname = paste0("./toyexample_nopr_sgd_pcn_nreal_", kreal, "_v6.npy")
            npySave(fname, nopr)
        }
    }
}




















