# PLScure
PLScure (which stands for <ins>**Partly Linear Single-index Cure**</ins> models with a
nonparametric, monotone increasing incidence link function) is a package that performs semiparametric estimation and inference for right-censored data using the method proposed by Lee et al. (2024).

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("https://github.com/lcyjames/PLScure/blob/main/PLScure.R?raw=TRUE")

# Usage #
The package contains 2 functions:
|Functions  | Description|
|------------- | -------------|
PLScureSIM  | Generate a data set according to the simulation study in Lee et al. (2024)
PLScureEST  | Perform the semiparametric estimation methods of Lee et al. (2024)

<ins>**PLScureSIM**</ins>
```
PLScureSIM(seed=NA, n, alpha, beta, gamma, scen)
```
This function generates a data set according to the model of the simulation study in Lee et al. (2024) that takes the arguments:
>- `n` is the sample size
>- `alpha` is the regression coefficient of `X` in the incidence component. If the Euclidean norm of this vector is not equal to 1, then the vector will be rescaled to norm 1.
>- `beta` is the regression coefficient of `W` in the latency component
>- `gamma` is the regression coefficient of `Z` in the latency component. If the Euclidean norm of this vector is not equal to 1, then the vector will be rescaled to norm 1.
>- `scen` is the setting used in the scenarios in the simulation study, which takes the values 1, 2, or 3 referring to Scenarios I, II, and III.

Example:
```
#This is the setting in Scenario III
Data <- PLScureSIM(seed = 1234, n = 500, alpha = c(1,-1,1,-1), beta = c(0.5,-0.5), gamma = c(1,-1), scen = 3)
head(Data)

#   id        Yi cen X1         X2         X3         X4 cure
# 1  1 1.6608376   0 -1  0.3115255  0.3143686  0.3592891    1
# 2  2 1.5471023   1  1  0.5060559 -0.5747400 -0.5466319    0
# 3  3 0.1586276   0 -1 -0.4771927 -0.9983864 -0.7762539    0
# 4  4 1.0180577   1 -1 -0.5110095 -0.9111954 -0.8371717    0
# 5  5 2.6545979   0  1 -0.4906859 -0.4405479  0.4595894    0
# 6  6 0.1450663   0 -1  0.5747557 -1.0236557 -0.0151383    1
```

This data structure is as follows:
>- `id` is the sample identifier
>- `Yi` is the exact failure time or censoring time
>- `cen` is the right-censoring indicator
>- `X1`, `X2`, `X3`, and `X4` are the covariates, where `X1` takes up the values -1 and 1 with equal probability, and `X2`, `X3` and `X4` are independent standard normal variables
>- `X1` and `X2` are included in covariate matrices `X` and `W` in the simulation setting
>- `X3` and `X4` are included in covariate matrices `X` and `Z` in the simulation setting
>- `cure` is the cure indicator, which is unobservable in real data analysis, hence not an input in the estimation function

<ins>**PLScureEST**</ins>

```
PLScureEST(X, W, Z, Yi, cen, K1 = 3, K2 = 5, M2 = 1:5, tolerance = 10^{-4}, attempt = 10, SE_est = TRUE, TRACE = TRUE)
```
This function performs the semiparametric estimation methods of Lee et al. (2024). The details of the arguments are as follows:
>- `X` is an n times p dimensional covariate matrix included in the single index of the incidence component; each column vector will be standardized to mean 0 and variance 1
>- `W` is an n times q dimensional covariate matrix included in the latency component with a linear effect on the log hazards
>- `Z` is an n times r dimensional covariate matrix included in the single index of the latency component; each column vector will be standardized to mean 0 and variance 1
>- `Yi` is a vector of observed event or censoring times
>- `cen` is a vector of censoring indicator that takes 0 for censoring and 1 for events
>- `K1` is the positive integer associated with m<sub>1</sub>, the number of basis functions used to approximate the `G` function, set to 3 by default
>- `K2` is the positive integer associated with the perturbation constant in the profile likelihood approach for the estimation of the standard error of Euclidean parameter estimates, set to 5 by default
>- `M2` is a vector of positive integers pertaining to the candidates for the search of the optimal value of m<sub>2</sub>, the number of basis functions used to approximate the `H` function, using AIC, set to 1:5 by default
>- `tolerance` is the stopping criterion for the EM algorithm, set to 10^{-4} by default
>- `attempt` is the number of random initializations of the parameter values used to perform the EM algorithm per each combination of `K1` and `M2`, set to 10 by default
>- `SE_est` is an option for computing the estimated standard error, set to TRUE by default
>- `TRACE` is an option for tracking the progress of the parameter estimation, set to TRUE by default

Example:
```
Data <- PLScureSIM(seed = 1234, n = 500, alpha = c(1,-1,1,-1), beta = c(0.5,-0.5), gamma = c(1,-1), scen = 3)
Result <- PLScureEST(X = cbind(Data$X1,Data$X2,Data$X3,Data$X4), W = cbind(Data$X1,Data$X2), Z = cbind(Data$X3,Data$X4),
                     Yi = Data$Yi, cen = Data$cen, K1 = 3, K2 = 5, M2 = 1:5, TRACE = TRUE)
Result

# $summary
#     coef estimate    SE CI.left CI.right
# 1 alpha1    0.551 0.176   0.206    0.897
# 2 alpha2   -0.498 0.159  -0.811   -0.186
# 3 alpha3    0.598 0.138   0.329    0.868
# 4 alpha4   -0.299 0.150  -0.593   -0.005
# 5  beta1    0.468 0.089   0.293    0.643
# 6  beta2   -0.537 0.089  -0.712   -0.363
# 7 gamma1    0.614 0.076   0.466    0.762
# 8 gamma2   -0.789 0.059  -0.905   -0.674
# 
# $alpha
# [1] 0.5512264 -0.4984034  0.5984616 -0.2993113
# 
# $alpha.se
# [1] 0.1762952 0.1592832 0.1375977 0.1500957
# 
# $beta
# [1] 0.4677577 -0.5371108
# 
# $beta.se
# [1] 0.08937872 0.08905294
# 
# $gamma
# [1] 0.6138690 -0.7894079
# 
# $gamma.se
# [1] 0.07559206 0.05878281
# 
# $psi
# [1] -4.7566623 -4.7560575 -4.7549804  0.3814506  0.3815638  0.3815658  0.3815853  0.3816561  0.3817943  0.3817962  0.3817983  0.3820711 33.3748788
# 
# $phi
# [1] 4.158311   7.133752 -23.451696  25.628132 -13.133059   4.073873
# 
# $m10.star
# [1] 12
# 
# $m20.star
# [1] 5
# 
# $Xmax
# 100% 
# 3.76596 
# 
# $Zmax
# 100% 
# 3.44924 
# 
# $likelihood
# [1] -1379.494
# 
# $AIC
# [1] 2812.989
# 
# $AIC.by.m10.m20
#    m10 m20     AIC
# 1  12   1 2848.032
# 2  12   2 2838.325
# 3  12   3 2822.233
# 4  12   4 2813.373
# 5  12   5 2812.989
```

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y., Wong, K. Y., and Bandyopadhyay D. (2024). Partly-linear single-index cure models with a nonparametric incidence link function. Statistical Methods in Medical Research [online], DOI:  [10.1177/09622802241227960](https://doi.org/10.1177/09622802241227960).
