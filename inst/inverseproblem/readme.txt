1) inverseproblem_second_moments_pcn.R
script for calculation of the second moments of increments.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

2) inverseproblem_sgd_pcn.R
script that performs MAP estimate of parameters with high resolution of forward solver and long SGD chains.
The MAP estimate of such estimators is considered as a ground-truth. 
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

3) inverseproblem_sgd_pcn.R
script that performs MAP estimate of parameters.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

4) inverseproblem_single_value_pcn.R
Script that generates several realizations of the estimator for the given quantity of interest.
Forward models with high resolution are utilized.
The mean value over the realizations is considered as a ground truth.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

5) inverseproblem_single_value_pcn.R
Script that generates several realizations of the estimator for the given quantity of interest.
Reflection maximal coupling of pCN proposals is utilized.
The data is generated and saved as .rds and .npy files

6) inverseproblem_solver_test.R
Script that generates computes square of the L2 norm of the differrence of forward model output
at subserqquent levels.
The data is generated and saved as .rds and .npy files

7) inverseproblem_data.npy
File that contains the data - vector of observations.

