# Generate data for analysis of convergence rate of the forward model
rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)
set.seed(1)

# data generation:
xdata_ref = npyLoad('toyexample_data.npy')  # load the data

su_param = 4.0  # std of prior distribution
sm_param = 1.0  # std of likelihood function

nu = 2  # dimension of model paramter space: nu = dimension!!!
u0 <- array(rep(0, nu), dim = c(nu))
u0[1] = 1.0  # ground truth value
u0[2] = -3.0  # ground truth value


# compute the difference between observation vectors at subsequent levels for different number of levels
nlevel = 10
res <- array(rep(0, nlevel), dim = c(nlevel))

for(klevel in 1 : nlevel)
{
    xobs = observation_toyexample(u0, klevel) - observation_toyexample(u0, klevel - 1)
    res[klevel] = sum(xobs * xobs)
}
npySave("./toyexample_solver_test.npy", res)







