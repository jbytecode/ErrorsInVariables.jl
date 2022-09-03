[![codecov](https://codecov.io/gh/jbytecode/eive.jl/branch/main/graph/badge.svg?token=KMF7H1DS01)](https://codecov.io/gh/jbytecode/eive.jl)

# Eives.jl
Error-in-variables estimation using Compact Genetic Algorithms in Julia

# The Problem 

Suppose the linear regression model is 

$$
y = \beta_0 + \beta_1 x^* + \varepsilon
$$

where $y$ is n-vector of the response variable, $\beta_0$ and $\beta_1$ are unknown regression parameteres, $\varepsilon$ is the iid. error term, $x^*$ is the unknown n-vector of the independent variable, and $n$ is the number of observations.

We call $x^*$ unknown because in some situations the true values of the variable cannot be visible or directly observable, or observable with some measurement error. Now suppose that $x$ is the observable version of the true values and it is defined as 

$$
x = x^* + \delta
$$

where $\delta$ is the measurement error and $x$ is the erroneous version of the true $x^*$. If the estimated model is 

$$
\hat{y} = \hat{\beta_0} + \hat{\beta_1}x 
$$

then the ordinary least squares (OLS) estimates are no longer unbiased and even consistent. 

Eive-cga is an estimator devised for this problem. The aim is to reduce the errors-in-variable bias with some cost of increasing the variance. At the end, the estimator obtains lower Mean Square Error (MSE) values defined as

$$
MSE(\hat{\beta_1}) = Var(\hat{\beta_1}) + Bias^2(\hat{\beta_1})
$$

for the Eive-cga estimator. For more detailed comparisons, see the original paper given in the Citation part. 

# Usage 

For the single variable case 

```Julia 
julia> eive(dirtyx = dirtyx, y = y, otherx = nothing) 
```

and for the multiple regression 

```Julia 
julia> eive(dirtyx = dirtyx, y = y, otherx = matrixofotherx) 
```

Note that the method assumes there is only one erroneous variable in the set of independent variables.

# Citation

```latex
@article{satman2015reducing,
  title={Reducing errors-in-variables bias in linear regression using compact genetic algorithms},
  author={Satman, M Hakan and Diyarbakirlioglu, Erkin},
  journal={Journal of Statistical Computation and Simulation},
  volume={85},
  number={16},
  pages={3216--3235},
  year={2015},
  publisher={Taylor \& Francis}
}
```