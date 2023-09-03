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
>- `alpha` is the regression coefficient of X in the incidence component. If the Euclidean norm of this vector is not equal to 1, then the vector will be rescaled to norm 1.
>- `beta` is the regression coefficient of W
>- `gamma` is the regression coefficient of Z in the latency component. If the Euclidean norm of this vector is not equal to 1, then the vector will be rescaled to norm 1.
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
>- `cure` is the cure indicator, which is unobservable in real data analysis, hence not an input in the estimation function

<ins>**PLScureEST**</ins>

```
PLScureEST(X, W, Z, Yi, cen, C1=3, C2=3, M2=1:5, tol=10^{-4}, attempt=10, SE_est=TRUE)
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
Data <- PLScureSIM(seed = 1234, n = 500, alpha = c(1,-1,1), beta = 0.5, gamma = c(1,-1), scen=2)
Result <- PLScureEST(X = cbind(Data$X1,Data$X2,Data$X3), W=Data$X1, Z=cbind(Data$X2,Data$X3), Yi = Data$Yi, cen = Data$cen, C1 = 3, C2 = 3, M2 = 1:5, tol = 10^{-4}, attempt = 10, SE_est = TRUE)
Result

# $summary
# coef estimate    SE CI.left CI.right
# 1 alpha1    0.622 0.051   0.522    0.721
# 2 alpha2   -0.654 0.062  -0.774   -0.533
# 3 alpha3    0.431 0.067   0.300    0.563
# 4  beta1    0.448 0.079   0.294    0.602
# 5 gamma1    0.773 0.050   0.674    0.871
# 6 gamma2   -0.635 0.061  -0.755   -0.515
# 
# $alpha
# [1]  0.6219258 -0.6535611  0.4313539
# 
# $alpha.se
# [1] 0.05076201 0.06162275 0.06698572
# 
# $beta
# [1] 0.4479409
# 
# $beta.se
# [1] 0.07863164
# 
# $gamma
# [1]  0.7725484 -0.6349559
# 
# $gamma.se
# [1] 0.05025599 0.06114628
# 
# $psi
# [1]  -4.2359447  -4.2156966  -0.2848389  -0.2839586  -0.2837657  -0.2836702  -0.2835927  -0.2834988  -0.2833435  -0.2829049
# [11]  -0.2811348  -0.2212790 564.8488788
# 
# $phi
# [1]  2.766668  4.839309 -7.565144  5.079797  2.947774
# 
# $m10.star
# [1] 12
# 
# $m20.star
# [1] 4
# 
# $Xmax
# 100% 
# 3.660898 
# 
# $Zmax
# 100% 
# 3.387188 
# 
# $likelihood
# [1] -1507.07
# 
# $AIC
# [1] 3062.14
# 
# $AIC.by.m10.m20
# m10 m20      AIC
# 1  12   1 3132.882
# 2  12   2 3094.603
# 3  12   3 3062.007
# 4  12   4 3058.140
# 5  12   5 3061.094
```

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y., Wong, K. Y., &Bandyopadhyay D. (2023+). Partly-linear single-index cure models with a nonparametric incidence link function. (Under review)
