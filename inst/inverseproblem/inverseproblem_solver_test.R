# This script generates data for estimation of the convergence rate of the forward solver
rm(list=ls())
library(Rcpp)
library(RcppCNPy)
library(UnbiasedMultilevel)

dimension = 2 # dimension of the system
xdata_ref = npyLoad("./inverseproblem_data.npy")  # load the generated data 
theta_param = 0.1  # theta parameter from the paper (1 / std^2) std - standard deviation in the likelihood

nu = 30 # number of model parameters for test
u0 <- array(rep(0, nu * dimension), dim = c(nu, dimension))
for(k in 1 : nu)
{
    u0[k, 1 : dimension] <- runif(dimension, min=-1.0, max=1.0)
}

nlevel = 10
res <- array(rep(0, nu * nlevel), dim = c(nu, nlevel))

for(k in 1 : nu)
{
    for(klevel in 1 : nlevel)
    {
        xobs = observation_inverseproblem(u0[k, 1 : dimension], klevel) - observation_inverseproblem(u0[k, 1 : dimension], klevel - 1)
        res[k, klevel] = sum(xobs * xobs)
    }
}
npySave("./inverseproblem_solver_test.npy", res)



