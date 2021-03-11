1) covid19_second_moments_pcn.R
script for calculation of the second moments of increments.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

2) covid19_sgd_pcn.R
script that performs MAP estimate of parameters with high resolution of forward solver and long SGD chains.
The MAP estimate of such estimators is considered as a ground-truth. 
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

3) covid19_sgd_pcn.R
script that performs MAP estimate of parameters.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

4) covid19_single_value_pcn.R
Script that generates several realizations of the estimator for the given quantity of interest.
Forward models with high resolution are utilized.
The mean value over the realizations is considered as a ground truth.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

5) covid19_single_value_pcn.R
Script that generates several realizations of the estimator for the given quantity of interest.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

6) covid19_solver_test.R
Script that computes square of the L2 norm of the differrence of forward model output
at subserqquent levels.
The data is generated and saved as .rds and .npy files

7) covid19_data.npy
File that contains the data - vector of observations.

8) covid19_integration.R
Script that computes expected value ov model predictions over model parameters for optimal (MAP) values of gamma
distribution parameters.

