# This script generates data to test that the variance of the estimator decreases as 1/N with number of realizations
rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)
set.seed(1)

xdata_ref = npyLoad("./covid19_data.npy")

dimension <- 3  # dimension of model parameters space
u0 <- array(rep(0, dimension), dim = c(dimension))
scaling_factor <- array(rep(1, dimension), dim = c(dimension))
# middle values for model parameters
u0[1] = 0.002
u0[2] = 0.300
u0[3] = 15.00
# scaling factor for model parameters
scaling_factor[1] = 0.001
scaling_factor[2] = 0.100
scaling_factor[3] = 10.00
umid_param = 0.0 + u0  
sf_param = 0.0 + scaling_factor

alpha_param = 1.0  # shape of the gamma distribution
beta_param = 1.0  # scale of the gamma distribution

log_lh_fun <- function(u, l, alpha_m = alpha_param, beta_m = beta_param,
                       xdata = xdata_ref, umid = umid_param, sf = sf_param)
{
    # output: log(likelihood) + log(prior)
    # u - vector of model parameters
    # l - level
    # alpha_m = alpha_param - shape of the gamma distribution
    # beta_m = beta_param - scale of the gamma distribution
    # umid = umid_param - shift of model parametsr
    # sf = sf_param - shift of model parameters
    # xdata = xdata_ref - vector of observations

    nu = length(u)  # dimesnion of model parameters
    umin = umid - sf  # calculation of the lower boundaries of model parameters
    umax = umid + sf  # calculation of the upper boundaries of model parameters

    ustar = 0.0 + umid  # valud of actual parameters of the ODE
    for(k0 in 1 : nu)
    {
        ustar[k0] = umid[k0] + sf[k0] * u[k0]
    }

    llh = 0.0
    # contribution of the prior diatribution
    for(k in 1 : length(umid))
    {
        if(umin[k] < ustar[k] && ustar[k] < umax[k] && is.finite(llh))
        {
            llh = llh + log(2.0 * sf[k])
        }
        else
        {
            llh = - Inf
        }
    }
    # contribution of the likelihood function
    if(is.finite(llh))
    {
        nmin = 29  # first day of observations
        nmax = 52  # last day of observations
        xmeas = observation_covid19(ustar, l)  # model predictions
        xmetric = log(xmeas[nmin : nmax]) - log(xdata[nmin : nmax])  # difference from true observations
        # log density of the gamma distribution
        for (k in 1 : length(xmetric))
        {
            if(xmetric[k] > 0 && is.finite(llh))
            {
                llh = llh + dgamma(xmetric[k], alpha_m, scale=beta_m, log=TRUE)
            }
            else
            {
                llh = - Inf
            }
        }
    }
    return(llh)
}


grad_log_lh_fun <- function(u, l, alpha_m = alpha_param, beta_m = beta_param,
                            xdata = xdata_ref, umid = umid_param, sf = sf_param)
{

    # output: d(log(likelihood) + log(prior)) / du
    # u - vector of model parameters
    # l - level
    # alpha_m = alpha_param - shape of the gamma distribution
    # beta_m = beta_param - scale of the gamma distribution
    # umid = umid_param - shift of model parametsr
    # sf = sf_param - shift of model parameters
    # xdata = xdata_ref - vector of observations

    nu = length(u)
    ustar = 0.0 + umid
    for(k in 1 : nu)
    {
        ustar[k] = ustar[k] + sf[k] * u[k]
    }
    xmeas = observation_covid19(ustar, l)
    xmeas_grad = observation_grad_covid19(ustar, l)
    grad_llh = 0.0 * u

    nmin = 29
    nmax = 52

    xmetric = log(xmeas[nmin : nmax]) - log(xdata[nmin : nmax])
    xmetric_grad = xmeas_grad[nmin : nmax, 1 : nu]
    for(k0 in 1 : nu)
    {
        for(k1 in 1 : length(xmetric))
        {
            xmetric_grad[k1, k0] = xmetric_grad[k1, k0] * sf[k0] / xmeas[nmin + k1 - 1]
        }
    }

    for (k0 in 1 : nu)
    {
        for (k in 1 : length(xmetric))
        {
            if(xmetric[k] > 0)
            {
                grad_llh[k0] = grad_llh[k0] + ((alpha_m - 1.0) / xmetric[k] - 1.0 / beta_m) * xmetric_grad[k, k0]
            }
        }
    }
    return(grad_llh)
}

objective_function <- function(u, l, alpha_m = alpha_param, beta_m = beta_param,
                               xdata = xdata_ref, umid = umid_param, sf = sf_param)
{
    # output: d(log(likelihood) + log(prior)) / dalpha_param,  d(log(likelihood) + log(prior)) / dbeta_param
    # u - vector of model parameters
    # l - level
    # alpha_m = alpha_param - shape of the gamma distribution
    # beta_m = beta_param - scale of the gamma distribution
    # umid = umid_param - shift of model parametsr
    # sf = sf_param - shift of model parameters
    # xdata = xdata_ref - vector of observations

    nu = length(u)
    nv = 2
    llh_grad <- array(rep(0, nv), dim = c(nv))  # memory allocation for the output variable
    ustar = 0.0 + umid  # calculations of ODE parameters
    for(k in 1 : nu)
    {
        ustar[k] = umid[k] + sf[k] * u[k]
    }

    nmin = 29
    nmax = 52
    xmeas = observation_covid19(ustar, l)  # model predictions
    xmetric = log(xmeas[nmin : nmax]) - log(xdata[nmin : nmax]) # argument of the Gamma distribution

    # derivative with repect to shape parameter (alpha)
    llh_grad_term = 0.0
    for(k in 1 : length(xmetric))
    {
        if(xmetric[k] > 0)
        {
            llh_grad_term = llh_grad_term + log(xmetric[k])
        }
    }
    llh_grad_term = llh_grad_term - length(xmetric) * digamma(alpha_m) - length(xmetric) * log(beta_m)
    llh_grad[1] = llh_grad[1] + llh_grad_term

    # derivative with repect to scale parameter (beta)
    llh_grad_term = 0.0
    llh_grad_term = sum(xmetric) / beta_m / beta_m - length(xmetric) * alpha_m / beta_m
    llh_grad[2] = llh_grad[2] + llh_grad_term
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

rinit <- function(level)
{
  # sample roughly from prior until you get the finine likelihood
  for(k in 1 : 1)
  {
    chain_state <- runif(dimension, min = -0.9, max = 0.9)
    current_pdf <- logtarget(level, chain_state)
    if(is.finite(current_pdf) == FALSE)
    {
        k = 1
    }
  }
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

beta = (8.0 + 1.0) / 2.0
w = init_pl(beta)

k0 = 500  # k from the paper
m0 = 2 * k0  # m form the paper
niter = 10080  # number os SGD steps
nreal = 1  # number of realisations
nv = 2  # dimension of the paramet space (of the gradient in sgd)

gf = 0.2  # factor for the SGD step-size. SGD step size is gf / sgd_step
xval = array(rep(0, nv, dim = c(nv)))


for(kreal in 1 : nreal)
{
    # memory allocation
    sgd_vec <- array(rep(0, niter * nv), dim = c(niter, nv))
    xval = array(rep(0, nv), dim = c(nv))
    nopr <- array(rep(0, niter), dim = c(niter)) # total operations counter
    alpha_param = 1.0  # shape of the gamma distribution
    beta_param = 1.0  # scale of the gamma distribution

    for(kiter in 1 : niter)
    {
        # sgd_vec = 0.0 * sgd_vec
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
                                          sampling_factor = min(0.5, 1.0 / (2.0 ^ (4 * vlevel + 1))),
                                          max_iterations = Inf)
            xval <- xval + x_increm$uestimator / w[1 + vlevel]
            nopr[kiter] = x_increm$cost
        }
        sgd_vec[kiter, 1 : nv] = 0.0 + xval[1 : nv]
        fname = paste0("./covid19_single_value_pcn_nreal_", kreal, ".rds")
        saveRDS(sgd_vec, file = fname)
        fname = paste0("./covid19_single_value_pcn_nreal_", kreal, ".npy")
        npySave(fname, sgd_vec)

        fname = paste0("./covid19_nopr_sv_pcn_nreal_", kreal, ".rds")
        saveRDS(nopr, file = fname)
        fname = paste0("./covid19_nopr_sv_pcn_nreal_", kreal, ".npy")
        npySave(fname, nopr)

    }
}



