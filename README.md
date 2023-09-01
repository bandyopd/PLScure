# PLScure
PLScure (which stands for <ins>**Partly Linear Single-Index cure**</ins> models with a
nonparametric incidence link function) is a package that performs semiparametric estimation and inference for right-censored data using the method proposed by Lee et al. (2023+).

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("...?raw=TRUE")

# Usage #
The package contains 2 functions:
|Functions  | Description|
|------------- | -------------|
PLScureSIM  | Generate a data set according to the simulation study in Lee et al. (2023+)
PLScureEST  | Perform the semiparametric estimation methods of Lee et al. (2023+)

<ins>**PLScureSIM**</ins>
```
PLScureSIM(seed=NA, n, alpha, beta, gamma, scen)
```
This function generates a data set according to the model of the simulation study in Lee et al. (2023+) that takes the arguments:
>- `n` is the sample size
>- `alpha` is the regression coefficient of X in the incidence component. If the norm of this vector is not equal to 1, then the vector will be normalized to norm 1.
>- `beta` is the regression coefficient of W
>- `gamma` is the regression coefficient of Z in the latency component. If the norm of this vector is not equal to 1, then the vector will be normalized to norm 1.
>- `scen` is the setting used in the scenarios in the simulation study, which takes the values 1, 2, or 3 referring to Scenarios I, II, and III.

Example:
```
#This is the setting in Scenario II
Data <- PLScureSIM(seed = 1234, n = 500, alpha = c(1,-1,1), beta = 0.5, gamma = c(1,-1), scen=2)
head(Data)

#  id         Yi cen         X1         X2         X3 cure
#1  1 8.71373045   0 -1.2070657  0.2774292  1.0844412    1
#2  2 1.50105375   1  0.5060559 -0.5747400 -0.5466319    0
#3  3 0.13600808   1 -0.7315360 -0.5166697 -1.7507334    0
#4  4 1.12204120   1  0.9594941 -0.1102855 -0.5110095    0
#5  5 1.05116810   1 -0.6470190  0.8681809  0.3756356    0
#6  6 0.04545912   0  0.4595894 -0.6937202 -1.4482049    0
```

This data structure is as follows:
>- `id` is the sample identifier
>- `Yi` is the exact failure time or censoring time
>- `cen` is the right-censoring indicator
>- `X1`, `X2`, and `X3` are the covariates, where each of them follows a standard normal distribution
>- `X1` is included in X and W
>- `X2` and `X3` are included in X and Z

<ins>**PLScureEST**</ins>

```
PLScureEST(data, P, m=10, tolerance=10^{-3}, gamma0=NA, beta0=NA, alpha10=NA, alpha20=NA, mu0=NA, sigma0=NA, TRACE=FALSE)
```
This function performs the semiparametric estimation methods of Lee and Wong (2023+). The details of the arguments are as follows:
>- `data` is a data.frame object shown in the above, with columns `id`, `Ti`, `cen`, `X[1]`,...,`X[P]`, `Z`
>- `P` is the dimension of covariate X, which is also equal to the dimension of gamma0
>- `m` is the number of nodes used in the Gaussian quadrature rule for truncated normal distributions
>- `tolerance` is the stopping criterion for the EM algorithm, set to 10^{-3} by default
>- `gamma0` is a vector of constants of size `P` for the initial values of parameter gamma, set to be rep(0,P) by default (gamma0=NA)
>- `beta0` is a constant for the initial value of parameter beta, set to be 0 by default (beta0=NA)
>- `alpha10` is a constant for the initial value of parameter alpha1, set to be 0 by default (alpha10=NA)
>- `alpha20` is a constant for the initial value of parameter alpha2, set to be 0 by default (alpha20=NA)
>- `mu0` is a constant for the initial value of parameter mu, set to be median of `Z` in `data` by default (mu0=NA)
>- `sigma0` is a constant for the initial value of parameter sigma, set to be 2 by default (sigma0=NA)
>- `TRACE` is an option for tracking the converging path of the parameter estimation, set to be FALSE by default

Example:
```
Dataset<-PLScureSIM(seed = 1234, n = 500, gamma = 0.5, beta = -1, alpha1 = 2, alpha2 = 1.5, mu = 1.5, sigma = 0.5)
Result <-RCPsurvEST(data = Dataset, P = 1, gamma0 = 0.5, beta0 = -1, alpha10 = 2, alpha20 = 1.5, mu0 = 1.5, sigma0 = 0.5,TRACE = F)
Result

# $loglik
# [1] -2456.817
# 
# $gamma.hat
# [1] 0.4504628
# 
# $beta.hat
# [1] -0.788148
# 
# $alpha1.hat
# [1] 1.751743
# 
# $alpha2.hat
# [1] 1.418709
# 
# $mu.hat
# [1] 1.603879
# 
# $sigma.hat
# [1] 0.3963745
# 
# $gamma.hat.se
# [1] 0.06527598
# 
# $beta.hat.se
# [1] 0.205753
# 
# $alpha1.hat.se
# [1] 0.2864374
# 
# $alpha2.hat.se
# [1] 0.2920645
# 
# $mu.hat.se
# [1] 0.1028265
# 
# $sigma.hat.se
# [1] 0.08603625
# 
# $conv
# [1] 0
```

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y., Wong, K. Y., &Bandyopadhyay D. (2023+). Partly-linear single-index cure models with a nonparametric incidence link function. (Under review)
