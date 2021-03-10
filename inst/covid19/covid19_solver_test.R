# This script generates data for estimation of the convergence rate of the forward solver
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

nlevel = 5
res <- array(rep(0, nlevel), dim = c(nlevel))

ustar = 0.0 + umid_param  # valud of actual parameters of the ODE

for(klevel in 1 : nlevel)
{
    xobs = observation_covid19(ustar, klevel) - observation_covid19(ustar, klevel - 1)
    res[klevel] = sum(xobs * xobs)
}
npySave("./covid19_solver_test.npy", res)



