# This script generates data to test that the variance of the estimator decreases as 1/N with number of realizations
rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)

xdata_ref = npyLoad("./inverseproblem_data.npy")  # load the generated data 
theta_param = 0.1  # theta parameter from the paper (1 / std^2) std - standard deviation in the likelihood


log_lh_fun <- function(u, l, theta_m = theta_param, xdata = xdata_ref, log_theta_mu = 0.0, log_theta_sd = 1.0)
{

    # output: log(likelihood) + log(prior)
    # u - vector of model parameters
    # l - level
    # theta parameter from the paper (1 / std^2) std - standard deviation in the likelihood
    # xdata = xdata_ref - data, or vector of observations
    # log-noraml prior on theta is supposed
    # log_theta_mu - mean value for log-normal prior on theta
    # log_theta_sd - standard deviation value for log-normal prior on theta

    nu = length(u)  # dimension of u
    umin <- array(rep(0, nu), dim = c(nu))  # lower boundaries of prior
    umax <- array(rep(0, nu), dim = c(nu))  # upper boundaries of prior
    for(k in 1 : nu)
    {
        umin[k] = -1.0
        umax[k] = 1.0
    }
    llh = 0.0  #  output variable
    out_of_domain = FALSE  # check if u is in the domain of prior distribution or not
    for (k in 1 : nu)
    {
        # if inside -> add contribution of the prior
        if((umin[k] < u[k]) && (u[k] < umax[k]))
        {
            llh = llh + log(umax[k] - umin[k])
        }
        # otherwise - set flag for out_of_domain
        else
        {
            out_of_domain = TRUE
        }
    }
    # if out_of_fomain - probability is zeros -> output log(0) = -Inf
    if(out_of_domain)
    {
        llh = -Inf
    }
    # otherwise - compute the contribution of the likelihood
    else
    {
        xmeas = observation_inverseproblem(u, l)
        # contribution of the likelihood
        for (k in 1 : length(xmeas))
        {
            llh = llh - 0.5 * log(2.0 * 3.1415926) + 0.5 * log(theta_m)
            llh = llh - 0.5 * theta_m * (xmeas[k] - xdata[k]) * (xmeas[k] - xdata[k])
        }
        # contribution of the prior distribution on theta: lognormal prior on theta
        llh = llh - 0.5 * log(2.0 * 3.1415926) - log(log_theta_sd)
        llh = llh - log(theta_m)
        llh = llh - (log(theta_m) - log_theta_mu) * (log(theta_m) - log_theta_mu) / log_theta_sd / log_theta_sd / 2.0
    }
    return(llh)
}

grad_log_lh_fun <- function(u, l, theta_m = theta_param, xdata = xdata_ref, log_theta_mu = 0.0, log_theta_sd = 1.0)
{

    # output: derivatives of (log(likelihood) + log(prior))
    # with respect to model parameters: u

    # u - vector of model parameters
    # l - level
    # theta parameter from the paper (1 / std^2) std - standard deviation in the likelihood
    # xdata = xdata_ref - vector of observations
    # log-noraml prior on theta is supposed
    # log_theta_mu - mean value for log-normal prior on theta
    # log_theta_sd - standard deviation value for log-normal prior on theta

    nu = length(u)  # dimension of model parameter vector
    xmeas = observation_inverseproblem(u, l)  # observation for vector of model parameters u and level l
    xmeas_grad = observation_grad_inverseproblem(u, l)  # derivatives of the numerical solution with respect to u at u and level l
    grad_llh = 0.0 * u  # init output

    for(k0 in 1 : nu)
    {
        for(k in 1 : length(xmeas))
        {
            grad_llh[k0] = grad_llh[k0] - theta_m * (xmeas[k] - xdata[k]) * xmeas_grad[k, k0]
        }
    }
    return(grad_llh)
}

objective_function <- function(u, l, theta_m = theta_param, xdata = xdata_ref, log_theta_mu = 0.0, log_theta_sd = 1.0)
{

    # gradient of log(Bayesian_evidence_factor) with respect to theta
    # output: dlog(Bayesian_evidence_factor) / dtheta_param

    # u - vector of model parameters
    # l - level
    # sigma_u = su_param - std of the prior distribution
    # sigma_m = sm_param - std of the likelihood
    # xdata = xdata_ref - vector of observations
    # log-noraml prior on theta is supposed
    # log_theta_mu - mean value for log-normal prior on theta
    # log_theta_sd - standard deviation value for log-normal prior on theta

    nu = length(u)
    nv = 3
    llh_grad <- array(rep(0, nv), dim = c(nv))
    xmeas = observation_inverseproblem(u, l)
    # contribution of the likelihood
    for (k in 1 : length(xmeas))
    {
        llh_grad[1] = llh_grad[1] + 0.5 * (1.0 / theta_m - (xmeas[k] - xdata[k]) * (xmeas[k] - xdata[k]))
    }
    # contribution of the prior on theta
    llh_grad[1] = llh_grad[1] - (log(theta_m) - log_theta_mu) / log_theta_sd / log_theta_sd / theta_m
    llh_grad[1] = llh_grad[1] - 1.0 / theta_m
    llh_grad[2] = u[1]
    llh_grad[3] = u[2]
    return(llh_grad)
}

init_pl <- function(beta, lmax=20)
{
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
proposal_sd = 1.0
proposal_rho = 0.95
tuning <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)
tuning_fine <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)
tuning_coarse <- list(proposal_sd = proposal_sd, proposal_rho = proposal_rho)

rinit <- function(level){
  chain_state <- runif(dimension, min=-0.9, max=0.9)
  current_pdf <- logtarget(level, chain_state)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

beta = (4.0 + 1.0) / 2.0  # value of P_L(l) decay
w = init_pl(beta)  # values of P_L(l)

k0 = 100  # k from the paper
m0 = 10 * k0  # m form the paper
niter = 10080  # number os SGD steps
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
    theta_param = 1.0  # std of likelihood function

    for(kiter in 1 : niter)
    {
        # sgd_vec = 0.0 * sgd_vec
        xval = 0.0 * xval
        vlevel = get_l_sample(w)
        if(vlevel == 0)
        {
            x_expect = unbiased_expectation(level = 0 + vlevel,
                                            rinit = rinit,
                                            single_kernel = single_pcn_kernel,
                                            tuning = tuning,
                                            coupled_kernel = coupled2_pcn_kernel,
                                            proposal_coupling = reflectionmaximal2_pcn_coupling,
                                            h = function(l, x) objective_function(x, l),
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
                                          h = function(l, x) objective_function(x, l),
                                          k = k0, m = m0,
                                          sampling_factor = min(0.5, 1.0 / (2.0 ^ (2 * vlevel + 1))),
                                          max_iterations = Inf)
            xval <- xval + x_increm$uestimator / w[1 + vlevel]
            nopr[kiter] = x_increm$cost
        }
        sgd_vec[kiter, 1 : nv] = 0.0 + xval[1 : nv]
        fname = paste0("./inverseproblem_single_value_pcn_nreal_", kreal, ".rds")
        saveRDS(sgd_vec, file = fname)
        fname = paste0("./inverseproblem_single_value_pcn_nreal_", kreal, ".npy")
        npySave(fname, sgd_vec)
        fname = paste0("./inverseproblem_nopr_sv_pcn_nreal_", kreal, ".rds")
        saveRDS(nopr, file = fname)
        fname = paste0("./inverseproblem_nopr_sv_pcn_nreal_", kreal, ".npy")
        npySave(fname, nopr)

    }
}

